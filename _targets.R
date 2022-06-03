library(targets)
library(tarchetypes)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Configuration ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tar_option_set(
  packages = c('lubridate',
               'stringr',
               'raster',
               'terra',
               'jsonlite',
               'tibble',
               'fst',
               'sf',
               'magrittr',
               'ggplot2',
               'SHAPforxgboost',
               'Just.universal',
               'xgboost',
               'parallel',
               'data.table'),
  format = 'qs',
  workspace_on_error = TRUE,
  error = 'abridge')

tar_config_set(store = '/data-coco/Earth_obs_cleaning/targets')
Sys.setenv(RSTUDIO_PANDOC = '/usr/lib/rstudio-server/bin/pandoc')
intermediate.path = function(...)
   file.path('/data-coco/Earth_obs_cleaning/intermediate', ...)
download = function(from, to, ...)
    download.update.meta(from, "/data-coco/Earth_obs_cleaning/downloads", to, ...)
satellite_hdf_root = '/data-coco/Earth_obs_cleaning/earthdata'

n.workers = 22L

source('R/globals.R')
source('R/data.R')
source('R/functions.R')
source('R/xgboost_cv.R')
source('R/paper_functions.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Targets ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

process_years = 2000:2021
all_dates = seq(
    lubridate::make_date(min(process_years)),
    lubridate::make_date(max(process_years), 12, 31),
    by = 1)
example_date = as.Date("2010-07-01")
  # This should be a date for which the satellite data of interest
  # exists on all tiles.
stopifnot(example_date %in% all_dates)

region_values = list(region = aoiname)
sat_values = list(sat = sats)
buffers_km = c(10, 30, 90, 270)
pred_round_digits = 5
  # MCD19 AOD has 3 digits of precision, so this is a little more.
agg_level = 10
agg_thresh = 3
features = c("MCD19_AOD_470nm", "dayint", "AOD_Uncertainty",
             "Column_WV", "cosSZA", "cosVZA", "RelAZ", "Scattering_Angle", "Glint_Angle", "qa_best",
             do.call(paste0, expand.grid(
               c("pNonNAAOD", "Mean_AOD", "diff_AOD"),
               paste0(buffers_km, "km"))))

terra.rast.fmt = tar_format(
    read = function(path)
        terra::rast(path),
    write = function(object, path)
        terra::writeRaster(object, path, filetype = "GTiff"))

set1_targets = list(
  tar_target(aer_stn_path,
             download(
                "https://aeronet.gsfc.nasa.gov/aeronet_locations_v3.txt",
                "aeronet_stations.txt")),
  tar_target(aer_files_path,
            {p = download(
                 "https://aeronet.gsfc.nasa.gov/data_push/V3/AOD/AOD_Level20_All_Points_V3.tar.gz",
                 "aeronet_observations.tar.gz")
             assert(0 == unlink(intermediate.path("aeronet"),
                 recursive = T))
             assert(dir.create(intermediate.path("aeronet")))
             assert(0 == system2("tar", shQuote(c(
                 "--extract", "--file", p,
                 "--directory", intermediate.path("aeronet")))))
             intermediate.path("aeronet/AOD/AOD20/ALL_POINTS/")}),
  tar_target(aer_stations,
             fread(aer_stn_path, col.names = c('Site_Name', 'lon', 'lat', 'elevm')),
             format = 'fst_dt'),
  tar_target(vrt_path,
             prepare_vrt_directory(intermediate.path())),

  # Area of interest ####
  tar_map( # region mapping
    values = region_values,
    tar_target(buff,
               get_aoi_buffer(region)),
    tar_target(satellite_hdf_files, get.earthdata(
               satellite_hdf_root,
               product = "MCD19A2.006",
               satellites = "terra.and.aqua",
               tiles = satellite_aod_tiles[[region]],
               dates = all_dates)),
    tar_target(pred_grid, format = terra.rast.fmt, make_pred_grid(
               satellite_hdf_files[date == example_date, path])),
    tar_target(ground_obs, format = "fst_dt",
               get_ground_obs(process_years, pred_grid)),

    # Load AERONET data ####
    tar_target(aer,
               select_stations(aer_stations, buff, terra::crs(
                   terra::rast(satellite_hdf_files[date == example_date, path[1]])))),
    tar_target(aer_nospace,
               sf::st_drop_geometry(aer)),
    tar_target(aer_data,
               get_stn_data(aod_dir = aer_files_path, stations = aer_nospace)),
    tar_target(aer_filtered,
               filter_aer_bydate(aer_data, all_dates),
               format = 'fst_dt'),

    # Load MCD19A2 AOD ####
    tar_map( # sat mapping
      values = sat_values,
      tar_target(mcd19_vars, derive_mcd19_vars(aer_filtered,
                                               n.workers = n.workers,
                                               load_sat = sat,
                                               buffers_km = buffers_km,
                                               aer_stn = as.data.table(aer),
                                               satellite_hdf_files = satellite_hdf_files,
                                               agg_level = agg_level,
                                               agg_thresh = agg_thresh,
                                               vrt_path = vrt_path),
                 format = 'fst_dt'),

      # Model ####
      tar_target(traindata,
                 prepare_dt(mcd19_vars, date_range = all_dates)),
      tar_map(
        values = list(loss = c("l1", "l2")),
        tar_target(initial_cv,
                   initial_cv_dart(traindata,
                                   absolute = loss == "l1",
                                   y_var = "diff_AOD",
                                   features = features,
                                   stn_var = "Site_Name"))),
      tar_target(full_model, dart_full(traindata,
                                       y_var = "diff_AOD",
                                       features = features)),

      # Prediction ####
      tar_map(
         # Branch on the year, as well as a fake year "test" for which
         # only one day is used.
         values = list(pred_year = c(
             list("test"), as.list(process_years))),

         tar_target(pred_out, format = "fst_dt",
                   {if (Sys.getenv("OMP_NUM_THREADS") != "1")
                      # https://github.com/dmlc/xgboost/issues/2094
                        stop("The environment variable OMP_NUM_THREADS must be set to 1 before R starts to avoid a hang in `predict.xgb.Booster`.")
                    rbindlist(parallel::mclapply(mc.cores = n.workers,
                        (if (pred_year == "test")
                            example_date else
                            all_dates[data.table::year(all_dates) == pred_year]),
                        function(this_date)
                            with.temp.seed(list("predict", this_date, sat), run_preds(
                                full_model, features,
                                grid = pred_grid,
                                round_digits = pred_round_digits,
                                data = pred_inputs(
                                    features = features,
                                    buffers_km = buffers_km,
                                    satellite_hdf_files = satellite_hdf_files,
                                    vrt_path = vrt_path,
                                    load_sat = sat,
                                    this_date = this_date,
                                    agg_level = agg_level,
                                    agg_thresh = agg_thresh,
                                    aoi = buff,
                                    pred_bbox = NULL)))))}),

        # Compare predictions to ground observations
        tar_target(ground_comparison, satellite_vs_ground(pred_out, ground_obs)),

        # Map Predictions ####
        tar_target(preds_ggplot,
                   ggplot_orig_vs_adj(pred_out, pred_dates[1], viz_op = 3),
                   packages = c('ggplot2', 'cowplot', 'data.table', 'fst')),
        tar_target(preds_mapshot,
                   mapshot_orig_vs_adj(pred_out, pred_dates[1], viz_op = 3,
                                       use_jenks = TRUE, maxpixels = 2e6),
                   packages = c('mapview', 'raster', 'data.table', 'fst', 'rgeoda'),
                   format = 'file')
      )
    )
  )
)

# Combine all initial CV results into a list ####
set1u = unlist(set1_targets)
combined_target = tar_combine(combined_cv, set1u[grep('initial_cv_l2', names(set1u))],
                              command = list(!!!.x))

# Render CV report ####
report_targets = list(
  tar_target(cv_summary_tables,
             cv_summary(combined_cv)),
  tar_render(initial_cv_report,
             'R/initial_cv_report.Rmd')
)

# Render CONUS AOD results ###
paper_conus_targets = list(
  tar_render(paper_conus_html,
             'R/CONUS_AOD.Rmd', output_format = "html_document",
             packages = c('sf', 'patchwork')),
    tar_render(paper_conus_pdf, output_format = "pdf_document",
             'R/CONUS_AOD.Rmd',
             packages = c('sf', 'patchwork'))
)
# Final targets list ####
list(set1_targets, combined_target, report_targets, paper_conus_targets)
