# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# * Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source('code/globals.R')
source('code/util.R')
source('code/data.R')
source('code/modeling.R')
source('code/paper_functions.R')

library(targets)

tar_option_set(
    packages = sapply(parse("code/libraries.R")[[1]][-1][[1]][-1],
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
if (!is.null(Wf$test_small_daterange) && Wf$test_small_daterange)
   {process_years = year(example_date)
    all_dates = example_date + (-1:1)}
stopifnot(example_date %in% all_dates)

terra.rast.fmt = tar_format(
  # None of `terra`'s output formats seems to round-trip properly,
  # so we implement our own.
    read = function(path)
       {o = qs::qread(path)
        r = terra::rast(
            nrows = o$rows,
            ncols = o$cols,
            ext = o$ext,
            crs = o$crs)
        for (lyr in names(o$values))
            suppressWarnings({r[[lyr]] = o$values[[lyr]]})
        r},
    write = function(object, path)
       {r = object
        qs::qsave(file = path, list(
            rows = nrow(r),
            cols = ncol(r),
            ext = as.vector(terra::ext(r)),
            crs = terra::crs(r),
            values = `names<-`(
                lapply(names(r), function(lyr)
                    as.data.frame(r[[lyr]], na.rm = F)[[1]]),
                names(r))))})

pbapply::pboptions(type = "timer")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# * Targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

list(

    # AERONET region-indepedent loading logic
    tar_target(aer_stn_path, download(
        "https://aeronet.gsfc.nasa.gov/aeronet_locations_v3.txt",
        "aeronet_stations.txt")),
    tar_target(aer_orig_obs_path, download(
        "https://aeronet.gsfc.nasa.gov/data_push/V3/AOD/AOD_Level20_All_Points_V3.tar.gz",
        "aeronet_observations.tar.gz")),
    tar_target(aer_files_path,
       {assert(0 == unlink(intermediate.path("aeronet"),
            recursive = T))
        assert(dir.create(intermediate.path("aeronet")))
        assert(0 == system2("tar", shQuote(c(
            "--extract", "--file", aer_orig_obs_path,
            "--directory", intermediate.path("aeronet")))))
        intermediate.path("aeronet/AOD/AOD20/ALL_POINTS/")}),
    tar_target(aer_stations, format = 'fst_dt', fread(
        aer_stn_path, col.names = c('Site_Name', 'lon', 'lat', 'elevm'))),

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
                tiles = satellite.tiles(buff),
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

    # AERONET processing for this region
    tar_target(aer, select_stations(
        aer_stations,
        buff,
        terra::crs(pred_grid))),
    tar_target(aer_filtered, format = 'fst_dt', filter_aer_bydate(
        dates = all_dates,
        get_stn_data(
            aod_dir = aer_files_path,
            stations = aer))),

    # This step is where most of the satellite data is read.
    tar_target(traindata, format = 'fst_dt', make_traindata(
        Wf$satellite.product, Wf$satellite,
        Wf$y.sat, Wf$features, Wf$window.radius,
        aer_filtered,
        aer,
        satellite_hdf_files,
        example_date,
        n.workers = pmin(8L, n.workers))),
          # Using a lot more workers on Coco seems to be slower.

    # Modeling
    tar_target(initial_cv, initial_cv_dart(
        traindata,
        y_var = "y.diff",
        features = Wf$features,
        stn_var = "site")),
    tar_target(full_model, dart_full(
        traindata,
        y_var = "y.diff",
        features = Wf$features)),

    # Summarize and report on the CV
    tar_target(cv_summary_table, cv.summary(
        initial_cv$mDT_wPred)),
    tarchetypes::tar_render(initial_cv_report,
        'writing/initial_cv_report.Rmd'),

    # Compare AQS to our predictions
    tar_target(aqs_obs, format = "fst_dt", get_aqs_obs(
        process_years, pred_grid)),
    tar_target(pred_at_aqs_sites, format = "fst_dt", new.preds.compact(
        dt.start = lubridate::as_datetime(min(aqs_obs$date) - 1),
        dt.end = lubridate::as_datetime(max(aqs_obs$date) + 1),
        cells = sort(unique(aqs_obs$cell)),
        targets = list(pred_grid, satellite_hdf_files, full_model))),
    tar_target(aqs_comparison, satellite_vs_aqs(
        pred_at_aqs_sites, aqs_obs)),

    # Data for maps
    tar_target(T.median.mse.map.data, median.mse.map.data(
        Wf$y.sat, Wf$satellite, Wf$satellite.product, config$n.workers,
        initial_cv$mDT_wPred, pred_grid, buff, satellite_hdf_files, full_model)),
    tar_target(T.baltimore.map.data, baltimore.map.data(
        pred_grid, buff, satellite_hdf_files, full_model)),

    # Render the CONUS AOD manuscript
    if (Wf$satellite.product == "mcd19a2") list(
        tarchetypes::tar_render(paper_conus_html, output_format = "html_document",
            packages = c('sf', 'patchwork'),
            'writing/CONUS_AOD.Rmd', quiet = F),
        tarchetypes::tar_render(paper_conus_pdf, output_format = "pdf_document",
            packages = c('sf', 'patchwork', 'magick'),
            'writing/CONUS_AOD.Rmd')))
