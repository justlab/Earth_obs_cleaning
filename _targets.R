library(targets)
library(tarchetypes)
library(future)
# library(future.callr)
# plan(callr)
# plan(multicore)
plan(multisession)

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
               'future',
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

# For speed, the following aren't dynamic files.
aer_files_path = '/data-coco/ECHO_PM/AeronetAODV3Level2/AOD/AOD20/ALL_POINTS/'
mcd19path = '/data-coco/mcd19/fst/conus_full'
hdf_root = '/mnt/qnap_geo/MCD19A2/HDF'

n.workers = 22L

source('R/globals.R')
source('R/data.R')
source('R/functions.R')
source('R/xgboost_cv.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Targets ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

process_years = 2000:2020
all_dates = seq(
    lubridate::make_date(min(process_years)),
    lubridate::make_date(max(process_years), 12, 31),
    by = 1)
region_values = list(regions = aoiname)
sat_values = list(sat = sats)
buffers_km = c(10, 30, 90, 270)
agg_level = 10
agg_thresh = 3
features = c("MCD19_AOD_470nm", "dayint", "AOD_Uncertainty",
             "Column_WV", "RelAZ", "qa_best",
             do.call(paste0, expand.grid(
               c("pNonNAAOD", "Mean_AOD", "diff_AOD"),
               paste0(buffers_km, "km"))))

set1_targets = list(
  tar_target(aer_stn_path,
             '/data-coco/ECHO_PM/AeronetAODV3Level2/AOD/AOD20/aeronet_locations_v3.txt',
             format = 'file'),
  tar_target(refgrid_path,
             '/data-belle/LST/MYD21A1/derived/conus_grid_2020.fst',
             format = 'file'),
  tar_target(refras_path,
             '/data-belle/LST/MYD21A1/derived/conus_myd21_stack.tif',
             format = 'file'),
  tar_target(aer_stations,
             fread(aer_stn_path, col.names = c('Site_Name', 'lon', 'lat', 'elevm')),
             format = 'fst_dt'),
  tar_target(vrt_path,
             prepare_vrt_directory(intermediate.path())),

  # Area of interest ####
  tar_map( # region mapping
    values = region_values,
    tar_target(buff,
               get_aoi_buffer(regions)),
    tar_target(aer,
               select_stations(aer_stations, buff, refgrid_path, refras_path)),
    tar_target(mcd_refras,
               crop_refras_mcd(refgrid_path, mcd19path, aoiname = regions)),

    # Load AERONET data ####
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
                                               hdf_root = hdf_root,
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
      tar_target(pred_dates, c(as.Date('2003-03-16'), as.Date('2003-03-17'))),

      tar_target(predinput, format = "fst_dt",
                {future::plan("multicore", workers = n.workers)
                 rbindlist(future.apply::future_lapply(pred_dates,
                     future.seed = c(terra = 1337, aqua = 1338)[sat],
                     function(this_date) pred_inputs(
                          features = features,
                          buffers_km = buffers_km,
                          hdf_root = hdf_root,
                          vrt_path = vrt_path,
                          load_sat = sat,
                          this_date = this_date,
                          agg_level = agg_level,
                          agg_thresh = agg_thresh,
                          aoi = buff,
                          pred_bbox = NULL)))}),

      tar_target(pred_out,
                 run_preds(predinput, full_model, features)),

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

# Combine all initial CV results into a list ####
set1u = unlist(set1_targets)
combined_target = tar_combine(combined_cv, set1u[grep('initial_cv', names(set1u))],
                              command = list(!!!.x))

# Render CV report ####
report_targets = list(
  tar_target(cv_summary_tables,
             cv_summary(combined_cv)),
  tar_render(initial_cv_report,
             'R/initial_cv_report.Rmd')
)

# Final targets list ####
list(set1_targets, combined_target, report_targets)
