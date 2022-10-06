# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# * Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source('R/globals.R')
source('R/data.R')
source('R/functions.R')
source('R/xgboost_cv.R')
source('R/paper_functions.R')

library(targets)
library(tarchetypes)

tar_config_set(
    config = file.path(data.dir, "targets.yaml"),
    store = file.path(data.dir, "targets_store"))
tar_option_set(
    packages = sapply(parse("R/libraries.R")[[1]][-1][[1]][-1],
        function(x) as.character(x[[2]])),
    format = 'qs',
    workspace_on_error = TRUE,
    error = 'abridge',
    memory = "transient",
    garbage_collection = TRUE)

daily.sat = Wf$satellite.product != "geonexl2"

example_date = switch(Wf$satellite.product,
  # This should be a date for which the satellite data of interest
  # exists on all tiles.
    mcd19a2 = as.Date("2010-07-03"),
    geonexl2 = as.Date("2018-07-03"))
process_years = switch(Wf$satellite.product,
    mcd19a2 = 2000 : 2021,
    geonexl2 = 2018 : 2019)

all_dates = seq(
    lubridate::make_date(min(process_years)),
    lubridate::make_date(max(process_years), 12, 31),
    by = 1)
if (Sys.getenv("EARTH_OBS_CLEANING_TEST_SMALL_DATERANGE") != "")
   {process_years = year(example_date)
    all_dates = example_date + (-1:1)}
stopifnot(example_date %in% all_dates)

buffers_km = c(10, 30, 90, 270)
pred_round_digits = 5
  # MCD19A2 AOD has 3 digits of precision, so this is a little more.
agg_level = 10
agg_thresh = 3
features = c(
    "MCD19_AOD_470nm", "dayint", "AOD_Uncertainty", "Column_WV",
    "cosSZA", "cosVZA", "RelAZ", "Scattering_Angle", "Glint_Angle",
    "qa_best",
    do.call(paste0, expand.grid(
        c("pNonNAAOD", "Mean_AOD", "diff_AOD"),
        paste0(buffers_km, "km"))))

terra.rast.fmt = tar_format(
    read = function(path)
        terra::rast(path),
    write = function(object, path)
        terra::writeRaster(object, path, filetype = "GTiff"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# * Targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set1_targets = list(

    # AERONET region-indepedent loading logic
    tar_target(aer_stn_path, download(
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
    tar_target(aer_stations, format = 'fst_dt', fread(
        aer_stn_path, col.names = c('Site_Name', 'lon', 'lat', 'elevm'))),

    tar_target(vrt_path, prepare_vrt_directory(
        intermediate.path())),

    tar_target(buff, get_aoi_buffer(
        Wf$region)),
    tar_target(satellite_hdf_files, switch(Wf$satellite.product,
        mcd19a2 = get.earthdata(
            satellite_hdf_root,
            product = "MCD19A2.006",
            satellites = (if (Wf$satellite %in% c("terra", "aqua"))
                "terra.and.aqua" else
                stop()),
            tiles = satellite_aod_tiles[[Wf$region]],
            dates = all_dates),
        geonexl2 =
           {paths = dir(file.path(geonexl2.dir, Wf$satellite),
                recursive = T, full.names = T)
            data.table(
                satellite = factor(Wf$satellite),
                datetime = as.POSIXct(strptime(tz = "UTC",
                    basename(paths),
                    "GO16_ABI12A_%Y%j%H%M_")),
                tile = factor(basename(dirname(paths))),
                path = paths)})),
    tar_target(pred_grid, format = terra.rast.fmt, make_pred_grid(
        Wf$satellite.product,
        satellite_hdf_files[
            (if (daily.sat)
                date == example_date else
                lubridate::as_date(datetime) == example_date),
            by = tile,
            head(.SD, 1)])),
    tar_target(ground_obs, format = "fst_dt", get_ground_obs(
        process_years, pred_grid)),

    # AERONET processing for this region
    tar_target(aer, select_stations(
        aer_stations,
        buff,
        terra::crs(pred_grid))),
    tar_target(aer_filtered, format = 'fst_dt', filter_aer_bydate(
        dates = all_dates,
        get_stn_data(
            aod_dir = aer_files_path,
            stations = sf::st_drop_geometry(aer)))),

    # This step is where most of the satellite data is read.
    tar_target(traindata, format = 'fst_dt', prepare_dt(
        date_range = all_date,
        derive_mcd19_vars(
            aer_filtered,
            n.workers = n.workers,
            load_sat = sat,
            buffers_km = buffers_km,
            aer_stn = as.data.table(aer),
            satellite_hdf_files = satellite_hdf_files,
            agg_level = agg_level,
            agg_thresh = agg_thresh,
            vrt_path = vrt_path))),

    # Modeling
    tar_map(values = list(loss = c("l1", "l2")),
        tar_target(initial_cv, initial_cv_dart(
            traindata,
            absolute = loss == "l1",
            y_var = "diff_AOD",
            features = features,
            stn_var = "Site_Name"))),
    tar_target(full_model, dart_full(
        traindata,
        y_var = "diff_AOD",
        features = features)),

    # Making new predictions
    tar_map(
        # Branch on the year, as well as a fake year "test" for which
        # only one day is used.
        values = list(pred_year = c(
            list("test"), as.list(process_years))),

        tar_map(values = list(pred_month = 1:12),
            tar_target(pred_out, format = "parquet",
               {if (Sys.getenv("OMP_NUM_THREADS") != "1")
                  # https://github.com/dmlc/xgboost/issues/2094
                    stop("The environment variable OMP_NUM_THREADS must be set to 1 before R starts to avoid a hang in `predict.xgb.Booster`.")
                rbindlist(parallel::mclapply(mc.cores = n.workers,
                    (if (pred_year == "test")
                        example_date else
                        all_dates[data.table::year(all_dates) == pred_year & data.table::month(all_dates) == pred_month]),
                    function(this_date) with.temp.seed(
                        list("predict", this_date, sat),
                        run_preds(
                            full_model,
                            features,
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
                                pred_bbox = NULL)))))}))),

    # Comparing AQS to our predictions
    tar_target(pred_at_aqs, format = "fst_dt", satellite_at_aqs_sites(
        Wf$region, process_years, sat, ground_obs)),
    tar_target(ground_comparison, satellite_vs_ground(
        pred_at_aqs, ground_obs)),

    # Maps
    tar_target(median_mse_date, initial_cv_l2$mDT_wPred
        [, by = aer_date, mean((diff_AOD_pred - diff_AOD)^2)]
        [which.min(abs(V1 - median(V1))), aer_date]),
    tar_target(preds_ggplot,
        packages = c('ggplot2', 'cowplot', 'data.table', 'fst', 'terra'),
        ggplot_orig_vs_adj(
            pred_out_4_2011[pred_date == as.Date("2011-04-18")],
            viz_op = 3,
            pred_grid)),
    tar_target(preds_mapshot, format = 'file',
        packages = c('mapview', 'raster', 'data.table', 'fst', 'rgeoda', 'terra'),
        mapshot_orig_vs_adj(
            pred_out_4_2011[pred_date == as.Date("2011-04-18")],
            viz_op = 3,
            pred_grid,
            use_jenks = TRUE,
            maxpixels = 2e6)))

# Render the CV report.
report_targets = list(
  tar_target(cv_summary_tables, cv_summary(
      initial_cv_l2)),
  tar_render(initial_cv_report,
      'R/initial_cv_report.Rmd'))

# Render CONUS AOD results.
paper_conus_targets = list(
    tar_render(paper_conus_html, output_format = "html_document",
        packages = c('sf', 'patchwork'),
        'R/CONUS_AOD.Rmd'),
    tar_render(paper_conus_pdf, output_format = "pdf_document",
        packages = c('sf', 'patchwork', 'magick'),
        'R/CONUS_AOD.Rmd'))

# The final targets list
list(set1_targets, report_targets, paper_conus_targets)
