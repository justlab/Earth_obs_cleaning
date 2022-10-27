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
    "y.sat", "time.sat",
    "y.sat.mean", "y.sat.present",
    "AOD_Uncertainty",
    "cosSZA", "cosVZA", "RelAZ", "Scattering_Angle", "Glint_Angle",
    (if (Wf$satellite.product == "mcd19a2")
        c("Column_WV", "qa_best")))

terra.rast.fmt = tar_format(
    read = function(path)
        terra::rast(path),
    write = function(object, path)
        terra::writeRaster(object, path, filetype = "GTiff"))

pbapply::pboptions(type = "timer")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# * Targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

list(

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
        mcd19a2 =
           {d = get.earthdata(
                satellite_hdf_root,
                product = "MCD19A2.006",
                satellites = (if (Wf$satellite %in% c("terra", "aqua"))
                    "terra.and.aqua" else
                    stop()),
                tiles = satellite_aod_tiles[[Wf$region]],
                dates = all_dates)
            setnames(d, "date", "time")
            d},
        geonexl2 =
          # As of 13 Oct 2022, the GeoNEX-L2 files are only available
          # at a location that's likely temporary, so proper automatic
          # download is not yet implemented.
           {paths = dir(file.path(geonexl2.dir, Wf$satellite),
                recursive = T, full.names = T)
            `[`(
                data.table(
                    satellite = factor(Wf$satellite),
                    time = as.POSIXct(tz = "UTC",
                        basename(paths),
                        "GO16_ABI12A_%Y%j%H%M_"),
                    tile = factor(basename(dirname(paths))),
                    path = paths),
                lubridate::as_date(time) %in% all_dates)})),
    tar_target(pred_grid, format = terra.rast.fmt, make_pred_grid(
        Wf$satellite.product,
        satellite_hdf_files[
            lubridate::as_date(time) == example_date,
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
    tar_target(traindata, format = 'fst_dt', make_traindata(
        Wf$satellite.product,
        Wf$satellite,
        features,
        aer_filtered,
        as.data.table(aer),
        satellite_hdf_files,
        example_date,
        n.workers = pmin(8L, n.workers))),
          # Using a lot more workers on Coco seems to be slower.

    # Modeling
    tar_map(values = list(loss = c("l1", "l2")),
        tar_target(initial_cv, initial_cv_dart(
            traindata,
            absolute = loss == "l1",
            y_var = "y.diff",
            features = features,
            stn_var = "site"))),
    tar_target(full_model, dart_full(
        traindata,
        y_var = "diff_AOD",
        features = features)),

    # Summarize and report on the CV
    tar_target(cv_summary_table, cv.summary(
        initial_cv_l2$mDT_wPred)),
    tar_render(initial_cv_report,
        'R/initial_cv_report.Rmd'))
