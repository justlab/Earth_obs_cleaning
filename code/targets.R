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

tar_file = tarchetypes::tar_file

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
            "--directory", intermediate.path("aeronet"),
            "--no-same-owner"))))
        intermediate.path("aeronet/AOD/AOD20/ALL_POINTS/")}),
    tar_target(aer_stations, format = 'fst_dt', fread(
        aer_stn_path, col.names = c('Site_Name', 'lon', 'lat', 'elevm'))),

    tar_target(region.shape, get.region.shape(
        Wf$region)),
    tar_target(satellite.files, switch(Wf$satellite.product,
        mcd19a2 =
           {d = get.earthdata(
                root.dir = satellite_hdf_root,
                bbox = st_bbox(st_transform(region.shape, crs.lonlat)),
                products = "MCD19A2_061",
                satellites = (if (Wf$satellite %in% c("terra", "aqua"))
                    "terra.and.aqua" else
                    stop()),
                tiles = satellite.tiles(region.shape),
                dates = Wf$dates)
            setnames(d, "date", "time")
            setkey(d[, -"product"], satellite, time, tile)[]},
        aodc = get.aodc(Wf$satellite, Wf$dates))),
    tar_target(pred.grid, format = terra.rast.fmt, get.pred.grid(
        Wf$satellite.product,
        region.shape,
        satellite.files[
            time == Wf$time.example,
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
        Wf$time.example,
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
    tarchetypes::tar_render(cv.report,
        'writing/cv_report.Rmd'),

    # Compare AQS to our predictions
    tar_target(aqs.obs, format = "fst_dt", get.aqs.obs(
        (if (daily.sat()) "day" else "hour"),
        Wf$years,
        pred.grid)),
    tar_target(pred.at.aqs.sites, format = "fst_dt", new.preds.compact(
        dt.start = (if (daily.sat())
            lubridate::as_datetime(min(aqs.obs$time) - 1) else
            min(aqs.obs$time)),
        dt.end = (if (daily.sat())
            lubridate::as_datetime(max(aqs.obs$time) + 1) else
            max(aqs.obs$time) + lubridate::hours(1)),
        cells = sort(unique(aqs.obs$cell)),
        targets = list(pred.grid, satellite.files, model.full))),
    tar_target(satellite.vs.aqs, get.satellite.vs.aqs(
        (if (daily.sat()) "day" else "hour"),
        pred.at.aqs.sites,
        aqs.obs)),
    if (Wf$satellite.product == "mcd19a2") list(
        tar_target(aqs.obs.yearly, format = "fst_dt", get.aqs.obs(
            "year", Wf$years, pred.grid)),
        tar_target(satellite.vs.aqs.yearly, get.satellite.vs.aqs(
            "year", pred.at.aqs.sites, aqs.obs.yearly))),

    # Data for maps
    tar_target(median.improve.map.data, get.median.improve.map.data(
        Wf$y.sat, Wf$satellite, Wf$satellite.product, config$n.workers,
        cv$mDT_wPred, pred.grid, region.shape, satellite.files, model.full)),
    tar_target(median.improve.map.center,
        get_conus()[as.integer(st_intersects(
            st_centroid(st_as_sfc(st_bbox(st_transform(crs = crs.lonlat,
                st_as_sf(
                     median.improve.map.data$pred[,
                         .(sat.coords.x, sat.coords.y)],
                     coords = c(1, 2),
                     crs = st_crs(tar_read(pred.grid))))))),
            st_transform(get_conus(), crs = crs.lonlat))),]),

    tar_target(baltimore.map.data, get.baltimore.map.data(
        pred.grid, region.shape, satellite.files, model.full)),
    tar_target(special.time.map.data,
       {time.min = lubridate::as_datetime("2023-06-07T12:00:00-0400")
        time.max = lubridate::as_datetime("2023-06-07T18:00:00-0400")
          # A time of intense smoke from Canadian wildfires in the US.
        corners = convert.crs(
            cbind(
                c(-81, -66),  # To the SW of Pennsylvania
                c( 39,  48)), # To the NE of CONUS
            crs.lonlat,
            terra::crs(pred.grid))
        cells = as.data.frame(terra::crop(tar_read(pred.grid), with(corners,
            terra::ext(x[1], x[2], y[1], y[2]))))$cell.local
        sat = satellite.files[time.min <= time & time <= time.max]
        ns = pbapply::pbsapply(cl = n.workers, sat$path, \(p)
            sum(!is.na(`[`(
                as.data.frame(terra::rast(p)$AOD, na.rm = F),
                cells, "AOD"))))
        dt = sat[which.max(ns), time]
        new.preds(
            dt, dt, cells = cells,
            targets = list(pred.grid, satellite.files, model.full))}),

    # Manuscript graphics
    tar_target(ipath,
       {idir = file.path(getwd(), "writing/img")
        dir.create(idir, showWarnings = F)
        \(name) file.path(idir, paste0(name, ".png"))}),
    tar_file(agreement.plot.path, ggsave(
        ipath("agreement_plot"), dpi = 300, width = 7, height = 4,
        agreement.plot(cv$mDT_wPred))),
    tar_target(shap_long, shap.prep(
        shap_contrib = cv$shap_score,
        X_train = cv$mDT_wPred[, mget(names(cv$shap_score))])),
    tar_file(shap.summary.plot.path, ggsave(
        ipath("shap_summary"), dpi = 300, width = 7, height = 7,
        shap.plot.summary(shap_long))),
    tar_file(shap.features.plot.path, ggsave(
        ipath("shap_features"), dpi = 300, width = 7, height = 7,
        shap.plot.dependence(data_long = tar_read(shap_long), x = 'time.sat', y = 'Column_WV', color_feature = 'Column_WV') +
            aes(x = lubridate::as_datetime(x_feature)) +
            ylab("SHAP value for MAIAC CWV") +
            theme(axis.title.x = element_blank(),
                  panel.grid.minor.x = element_blank()))),
    tar_file(pred.map.median.improve.path, ggsave(
        ipath("pred_map_median_improve"), width = 7, height = 7,
        pred.map(
            median.improve.map.data$pred, pred.grid,
            bg.sf = get_conus(), color.scale.name = "AOD",
            quantile.cap = .99))),
    tar_file(pred.map.baltimore.path, ggsave(
        ipath("pred_map_baltimore"), width = 7, height = 7,
        pred.map(
            baltimore.map.data, pred.grid,
            bg.sf = get_conus(), color.scale.name = "AOD"))),
    tar_file(pred.map.special.time.path, ggsave(
        ipath("pred_map_special_time"), width = 7, height = 7,
        pred.map(
            special.time.map.data, pred.grid,
            bg.sf = get_conus(), color.scale.name = "AOD"))),

    # Render the manuscript.
    tarchetypes::tar_quarto(paper,
        sprintf('writing/paper_%s.qmd', Wf$satellite.product), quiet = F),

    # Do a supplementary analysis for AODC.
    if (Wf$satellite.product == "aodc") list(
        tar_target(zhang.comparison.data, get.zhang.comparison.data(
            Wf$y.sat,
            aer_filtered, aer, satellite.files, pred.grid,
            keep.qualities = c(0L, 1L),
            n.workers = pmin(8L, n.workers))),
        tar_target(zhang.comparison, do.zhang.comparison(
            zhang.comparison.data))))
