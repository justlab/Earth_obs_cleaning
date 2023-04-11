## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## * Setup
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## * Targets
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

    tar_target(region.shape, get.region.shape(
        Wf$region)),
    tar_target(satellite.files, switch(Wf$satellite.product,
        mcd19a2 =
           {d = get.earthdata(
                satellite_hdf_root,
                product = "MCD19A2.006",
                satellites = (if (Wf$satellite %in% c("terra", "aqua"))
                    "terra.and.aqua" else
                    stop()),
                tiles = satellite.tiles(region.shape),
                dates = Wf$dates)
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
                lubridate::as_date(time) %in% Wf$dates)})),
    tar_target(pred.grid, format = terra.rast.fmt, get.pred.grid(
        Wf$satellite.product,
        satellite.files[
            lubridate::as_date(time) == Wf$date.example,
            by = tile,
            head(.SD, 1)])),

    # AERONET processing for this region
    tar_target(aer, select_stations(
        aer_stations,
        region.shape,
        terra::crs(pred.grid))),
    tar_target(aer_filtered, format = 'fst_dt',
        get_stn_data(
            aod_dir = aer_files_path,
            stations = aer)[aer_date %in% Wf$dates]),

    # This step is where most of the satellite data is read.
    tar_target(traindata, format = 'fst_dt', get.traindata(
        Wf$satellite.product, Wf$satellite,
        Wf$y.sat, Wf$features, Wf$window.radius,
        aer_filtered,
        aer,
        satellite.files,
        Wf$date.example,
        n.workers = pmin(8L, n.workers))),
          # Using a lot more workers on Coco seems to be slower.

    # Modeling
    tar_target(cv, cv_dart(
        traindata,
        y_var = "y.diff",
        features = Wf$features,
        stn_var = "site")),
    tar_target(model.full, dart_full(
        traindata,
        y_var = "y.diff",
        features = Wf$features)),

    # Summarize and report on the CV
    tar_target(cv.summary, get.cv.summary(
        cv$mDT_wPred)),
    tarchetypes::tar_render(cv_report,
        'writing/cv_report.Rmd'),

    # Compare AQS to our predictions
    tar_target(aqs.obs, format = "fst_dt", get.aqs.obs(
        Wf$years, pred.grid)),
    tar_target(pred.at.aqs.sites, format = "fst_dt", new.preds.compact(
        dt.start = lubridate::as_datetime(min(aqs.obs$date) - 1),
        dt.end = lubridate::as_datetime(max(aqs.obs$date) + 1),
        cells = sort(unique(aqs.obs$cell)),
        targets = list(pred.grid, satellite.files, model.full))),
    tar_target(satellite.vs.aqs, get.satellite.vs.aqs(
        pred.at.aqs.sites, aqs.obs)),

    # Data for maps
    tar_target(median.mse.map.data, get.median.mse.map.data(
        Wf$y.sat, Wf$satellite, Wf$satellite.product, config$n.workers,
        cv$mDT_wPred, pred.grid, region.shape, satellite.files, model.full)),
    tar_target(baltimore.map.data, get.baltimore.map.data(
        pred.grid, region.shape, satellite.files, model.full)),

    # Render the CONUS AOD manuscript
    if (Wf$satellite.product == "mcd19a2")
        tarchetypes::tar_quarto(paper_conus_html,
            'writing/CONUS_AOD.qmd', quiet = F))
