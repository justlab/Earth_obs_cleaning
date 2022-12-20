# Contains most functions required for the AOD cleaning targets workflow

#' Return AERONET points that intersect the polygon describing the area of
#' interest.
#'
#' @param stations data.table including coordinates (\code{lon, lat}) points of
#'   AERONET stations
#' @param reg_polygon SF polygon defining the area of interest
#' @param new_crs the CRS to set the result to
#' @return a subset of geometry points within the given polygon
select_stations <- function(stations, reg_polygon, new_crs){
  stations = stations[
      Site_Name != "CART_SITE" &
        # This station is in the same place as "Cart_Site" and has
        # relatively few observations. It's not clear to what degree it's
        # a duplicate.
      in.sf(lon, lat, crs.lonlat, reg_polygon)]
  setkey(stations, Site_Name)
  stations[, c("x", "y") :=
      convert.crs(cbind(lon, lat), crs.lonlat, new_crs)]
  stations
}

#' Load Aeronet AOD measurement data from text files sourced from downloaded
#' \code{tar.gz} file. Only load the data for the specified stations.
#'
#' @param aod_dir directory where the \code{.lev20} files were extracted to.
#' @param stations SF points of AERONET stations
#' @param date_start subset to observations this date or later
#' @param date_end subset to observations this date or earlier
#' @return data.table of AERONET observations, joined with MODIS reference grid
#'   unique ID
get_stn_data <- function(aod_dir, stations, date_start = NULL, date_end = NULL){
  assert(!anyDuplicated(stations$Site_Name))
  # open files containing names of stations in the specified region
  aer_files_dir = list.files(aod_dir, full.names = TRUE)
  aer_files_dir = aer_files_dir[
    str_match(basename(aer_files_dir),
      "\\A\\d{8}_\\d{8}_(.+)\\.lev20\\z")[,2]
    %in% stations$Site_Name]
  if(length(aer_files_dir) > 0){
    t0 <- fread(aer_files_dir[[1]], nrows = 10) # to get variable names:
    vars_aod <- intersect(grep("AOD_", names(t0), value = T), grep("nm", names(t0), value = T))
    # sort the wave lengths varnames from low to high, and update the vars_aod
    x_nm <- sort(readr::parse_number(vars_aod)) # 340, 380, ... , 1640
    vars_aod <- paste0("AOD_", x_nm,"nm") # the AOD_...nm variables
    vars_wv <- paste0("Exact_Wavelengths_of_AOD(um)_", x_nm,"nm") # the Exact wv length

    # data_path is a single file path
    read_aod <- function(data_path, date_start = NULL, date_end = NULL){
      dt <- fread(data_path, select =  c(vars0, vars_aod, vars_wv))
      dt[, stn_time := as.POSIXct(paste(`Date(dd:mm:yyyy)`, `Time(hh:mm:ss)`),
                                  format = "%d:%m:%Y %H:%M:%S", tz = "UTC")]
      dt[, c("Date(dd:mm:yyyy)", "Time(hh:mm:ss)") := NULL]
      for (i in names(dt)) dt[get(i) == -999, (i):= NA] # set NA
      if(!is.null(date_start)) {dt = dt[as.Date(stn_time) >= as.Date(date_start), ] }
      if(!is.null(date_end))   {dt = dt[as.Date(stn_time) <= as.Date(date_end), ] }
      setcolorder(dt, c("AERONET_Site_Name", "stn_time"))
      dt
    }

    file_list <- lapply(unlist(aer_files_dir), read_aod, date_start = date_start, date_end = date_end)
    aer_data <- rbindlist(file_list)

    aer_data[, aer_date := as.Date(stn_time)]
    assert(nrow(unique(aer_data[, .(AERONET_Site_Name, stn_time)])) ==
      nrow(aer_data))
    aer_data
  } else {
    NULL
  }
}

#' Filter data.table of AERONET observations to only those in the given dates
#'
#' @param aer_data AERONET observation data.table with a date column `aer_date`
#' @param dates vector of dates to filter the AERONET observation
#' @return data.table of AERONET observations in the given dates
#'
filter_aer_bydate <- function(aer_data, dates){
  aer_data[aer_date %in% dates, ]
}

# Interpolate AERONET AOD to the desired wavelength.
interpolate_aod <- function(aer_data)
   {# N.B. All computations with wavelengths here are in micrometers,
    # even though the column names of the observations data use
    # nanometers, as in "AOD_531nm". We do this because the input
    # exact wavelengths are in micrometers.
    target.wl = 0.47
    max.input.wl = 1.0
      # This value was chosen because AOD at wavelengths above 1 μm
      # seemed to lie off the trendline for interpolation of 0.47 μm.
    min.obs = 4

    # Reshape the observations to one row per (site, time, wavelength,
    # AOD).
    d = `[`(
        melt(aer_data, measure.vars = patterns(
            exact.wl = "Exact_Wavelength",
            aod = "AOD_")),
        exact.wl <= max.input.wl & !is.na(aod) & aod > 0,
        -"variable")
    assert(!anyNA(d))

    # Run the interpolation models.
    d[, keyby = .(site, time.ground),
        if (.N >= min.obs &
                # Ensure we have at least one observation at
                # wavelengths less than and greater than the target.
                any(exact.wl <= target.wl) &
                any(exact.wl >= target.wl))
           {m = lm(log(aod) ~ log(exact.wl) + I((log(exact.wl))^2))
            if (m$rank == length(coef(m)))
              # Only use the model if it's full-rank.
                list(aod = exp(predict(m,
                    data.frame(exact.wl = target.wl))))}]}

make_traindata = function(
        satellite.product, the.satellite,
        y.sat.name, features, window.radius,
        aer_filtered, aer_stn, satellite_hdf_files, example_date,
        n.workers)
   {temporal.matchup.seconds = (if (daily.sat(satellite.product))
        8 * 60 else
        30)

    message("Assembling ground data")
    aer_filtered = aer_filtered[, c(
        list(
            site = AERONET_Site_Name,
            time.ground = stn_time),
        mget(str_subset(colnames(aer_filtered), "AOD")))]
    setkey(aer_filtered, site, time.ground)
    ground = aer_filtered[, .(site, time.ground)]
    ground[, c("x_satcrs", "y_satcrs") := aer_stn[
        .(ground$site), on = "Site_Name", .(x, y)]]
    setkey(ground, time.ground)

    satellite_hdf_files = copy(satellite_hdf_files)
    satellite_hdf_files[, sat.files.ix := .I]

    all.tiles = unique(satellite_hdf_files$tile)
    d = rbindlist(lapply(all.tiles, function(the.tile)
       {message(sprintf("Getting training data for tile %s (%d / %d)",
            the.tile, match(the.tile, all.tiles), length(all.tiles)))

        # Set satellite cell indices for the ground data, dropping
        # ground data outside of this tile.
        r.example = read_satellite_raster(satellite.product, the.tile,
            satellite_hdf_files[j = path[1],
                tile == the.tile &
                lubridate::as_date(time) == example_date])
        d = cbind(ground, cell.local = terra::cellFromXY(
            r.example,
            ground[, .(x_satcrs, y_satcrs)]))[!is.na(cell.local)]
        if (!nrow(d))
            return()

        # Identify the satellite data for this tile.
        sat = expand.to.overpasses(
            satellite_hdf_files[tile == the.tile,
                .(time.sat = time, sat.files.ix)],
            satellite_hdf_files,
            the.satellite, satellite.product, n.workers)
        if (!nrow(sat))
            return()
        assert(!anyNA(sat[, -"overpass"]))
        setkey(sat, time.sat)

        # Join the satellite files with the ground data by time.
        message("Roll-joining")
        d = cbind(d[, -"time.ground"], sat[d,
            roll = "nearest", mult = "first",
            .(time.ground = i.time.ground, time.sat = x.time.sat,
                sat.files.ix, overpass)])
        d = d[!is.na(time.sat)]
        d[, time.diff := as.numeric(difftime(time.ground, time.sat,
            units = "secs"))]
        d = d[abs(time.diff) <= temporal.matchup.seconds]
        if (!nrow(d))
            return()

        # Use each satellite cell-time only once, taking the closest
        # time matchup.
        message("Reducing")
        d = d[order(abs(time.diff)), by = .(time.sat, cell.local),
            head(.SD, 1)]

        # Read in the values of the appropriate satellite files.
        get.predictors(
            d, satellite_hdf_files,
            satellite.product, y.sat.name, features, window.radius, n.workers)}))

    message("Interpolating AOD")
    d[, y.ground := `[`(
        interpolate_aod(aer_filtered[.(d$site, d$time.ground)]),
        .(d$site, d$time.ground),
        aod)]
    d = d[!is.na(y.ground)]

    # Calculate the difference.
    d[, y.diff := y.sat - y.ground]
    d}

expand.to.overpasses = function(d, satellite_hdf_files, the.satellite, satellite.product, n.workers)
   {if (!multipass.sat(satellite.product))
        return(cbind(d, overpass := NA_integer_))
    # Each file contains multiple overpasses and multiple
    # satellites, stored in different layers. Expand `d` into
    # one row per layer, keeping only the satellite of interest.
    message("Collecting overpasses")
    rbindlist(pblapply(cl = n.workers, 1 : nrow(d), function(i)
       {path = satellite_hdf_files[d[i, sat.files.ix], path]
        if (file.size(path) == 0)
            return()
        orbit.times = str_extract_all(
            fromJSON(terra::describe(options = "json", path))
                $metadata[[1]]$Orbit_time_stamp,
            "\\S+")[[1]]
        `[`(
            d[i, .(
                sat.files.ix,
                overpass = seq_along(orbit.times),
                satellite = c(T = "terra", A = "aqua")[
                    str_sub(orbit.times, -1)],
                time.sat = as.POSIXct(tz = "UTC",
                    orbit.times,
                    "%Y%j%H%M"))],
            satellite == the.satellite,
            -"satellite")}))}

get.predictors = function(
        d, satellite_hdf_files,
        satellite.product, y.sat.name, features, window.radius, n.workers)
   {vnames = intersect(
        str_replace(features, "qa_best", "AOD_QA"),
        feature.raster.layers)
    window.matrix = as.matrix(CJ(
        row = -window.radius : window.radius,
        col = -window.radius : window.radius))
    chunks = split(d, by = c("sat.files.ix", "overpass"))
    parallel.outside = length(chunks) > 5
    xlapply = function(do.pblapply, ...)
        if (do.pblapply)
            pblapply(cl = n.workers, ...)
        else
            lapply(...)

    message("Reading predictors from satellite data")
    d = rbindlist(xlapply(parallel.outside, chunks, \(chunk) chunk
        [, c("y.sat", vnames, "y.sat.mean", "y.sat.present") :=
           {r = read_satellite_raster(
                satellite.product,
                satellite_hdf_files[chunk$sat.files.ix[1], tile],
                satellite_hdf_files[chunk$sat.files.ix[1], path],
                overpass[1])
            y.sat.values = drop(r[[y.sat.name]][])
            cbind(
                y.sat.values[cell.local],
                r[[vnames]][cell.local],
                rbindlist(xlapply(!parallel.outside, cell.local, \(cell)
                   {rc = terra::rowColFromCell(r, cell)
                    v = y.sat.values[na.omit(cellFromRowCol(r,
                        rc[,1] + window.matrix[,1],
                        rc[,2] + window.matrix[,2]))]
                    data.frame(
                        y.sat.mean = mean(v, na.rm = T),
                        y.sat.present = mean(!is.na(v)))})))}]
        [!is.na(y.sat)]))
          # Keep only cases where we have the satellite
          # outcome of interest.

    if ("time.diff" %in% colnames(d))
        d[, time.diff := NULL]

    if ("AOD_QA" %in% colnames(d))
       {# See page 13 of https://web.archive.org/web/20200927141823/https://lpdaac.usgs.gov/documents/110/MCD19_User_Guide_V6.pdf
        # Some variables produced here may not be used in training,
        # but may still be used for stratifying CV results.
        bits = function(x, bit.index.lo, bit.index.hi)
            bitwAnd(
                bitwShiftL(1L, bit.index.hi - bit.index.lo + 1L) - 1L,
                bitwShiftR(x, bit.index.lo))
        d[, qa_land := structure(class = "factor",
            bits(AOD_QA, 3L, 4L) + 1L,
            levels = c("Land", "Water", "Snow", "Ice"))]
        d[, qa_adjacent := structure(class = "factor",
            bits(AOD_QA, 5L, 7L) + 1L,
            levels = c(
                "Normal", "AdjacentClouds", "Surrounded",
                "AdjacentSingleCloudy", "AdjacentSnow", "PreviousSnow"))]
        assert(d[, all(as.integer(qa_adjacent) <= 5L)])
        d[, qa_best := bits(AOD_QA, 8L, 11L) == 0L]
        d[, qa_model := structure(class = "factor",
            bits(AOD_QA, 13L, 14L) + 1L,
            levels = c("Background", "Smoke", "Dust"))]
        assert(d[, all(as.integer(qa_model) <= 3L)])
        d[, AOD_QA := NULL]}

    d}

#' Do initial CV with hyperparameter selection with DART and rank features.
#'
#' Must provide either stn_var or day_var for binning folds.
#'
#' @param data data.frame of observations with DV and IVs
#' @param y_var column name to predict
#' @param features character vector of column names to train on
#' @param stn_var if binning by stations, provide the station column name
#' @param day_var if binning by days, provide the day column name
# Depends on functions in xgboost_cv_RFE.R
initial_cv_dart <- function(
  data,
  y_var,
  features,
  k_fold = 5,
  n_rounds = 100,
  stn_var = NULL,
  day_var = NULL,
  progress = TRUE
){
  xgb_threads <- get.threads()

  mDT <- data.table::copy(setDT(data)) # needed for drake, uncertain about targets
  simplify.dt.for.xgboost(mDT)
  if(!is.null(day_var)) {
    mDT[, dayint := as.integer(get(day_var))]
    by_var = "day"
  }
  if(!is.null(stn_var)) {
    mDT[, stn := get(stn_var)]
    by_var = "stn"
  }
  if (!is.null(day_var) & !is.null(stn_var))
    stop("Please provide either day_var or stn_var.")

  # some internal variables names
  y_formula <- as.formula(paste(y_var, "~."))
  y_var_pred <- paste0(y_var, "_pred") # name of the predicted y
  y_var_pred_whole <- paste0(y_var, "_pred_whole") # name of the predicted y

  # bin data into specified number of folds
  temp <- prepare.bin(mDT, by_var = by_var, k_fold = k_fold)
  mDT <- temp$data
  bin_list <- temp$bin # list of values of selected variable (stn or day) by fold
  rm(temp)
  index.fs.list <-  list(
    # return the observations in each fold
    stn = function(x, bin_list) mDT[stn%in%bin_list[[x]], which = TRUE],
    day = function(x, bin_list) mDT[dayint%in%bin_list[[x]], which = TRUE]
  )

  index_train <- index_test <- list()
  for (k in 1:k_fold){
    index_test[[k]]  <- index.fs.list[[by_var]](k, bin_list)
    index_train[[k]] <- -index_test[[k]]
  }

  # run k-fold cv and record SHAP matrix, predicted value, overall rmse...
  if (!y_var%in%features) features <- c(features, y_var)

  message("Run k-fold cv \n")
  cv_results <- run.k.fold.cv(k_fold = k_fold,
                              dataXY_df = mDT[, ..features],
                              by_var = by_var,
                              n_rounds = n_rounds,
                              y_var = y_var,
                              progress = progress,
                              seed = 1234,
                              index_train = index_train,
                              index_test = index_test,
                              xgb_threads = xgb_threads)

  mDT <- cbind(mDT, cv_results$y_pred_dt)
  mDT[, time.sat := lubridate::as_datetime(time.sat)]
  mDT[, y.ground.pred := y.sat - y.diff_pred]

  # output
  c(list(by_var = by_var, bin_list = bin_list, mDT_wPred = mDT), cv_results)
}

#' Summarize CV statistics.
#'
#' @param d The data table `mDT_wPred` returned by `initial_cv_dart`.
cv.summary = function(d)
    rbind(
        d[, c(list(Year = "all"), eval(performance.j))],
        d[, keyby = .(Year = year(time.sat)), eval(performance.j)])

performance.j = quote(
   {mse = function(x, y) mean((x - y)^2)
    sdn = function(x) sqrt(mse(x, mean(x)))
    tqs = unname(quantile(
        abs(as.numeric(difftime(units = "mins",
            time.sat, time.ground))),
        c(.25, .5, .75)))
    list(
        "Cases" = .N,
        "Sites" = length(unique(site)),
        "RMSE, raw" = sqrt(mse(y.ground, y.sat)),
        "RMSE, corrected" = sqrt(mse(y.ground, y.ground.pred)),
        "Proportion of raw MSE" = mse(y.ground, y.ground.pred) / mse(y.ground, y.sat),
        "Median, ground" = median(y.ground),
        "SD, ground" = sdn(y.ground),
        "SD, raw" = sdn(y.sat),
        "SD, corrected" = sdn(y.ground.pred),
        "Bias, raw" = mean(y.sat - y.ground),
        "Bias, corrected" = mean(y.ground.pred - y.ground),
        "Correlation, raw" = cor(y.sat, y.ground),
        "Correlation, corrected" = cor(y.ground.pred, y.ground),
        "Time difference, Q1" = tqs[1],
        "Time difference, median" = tqs[2],
        "Time difference, Q3" = tqs[3])})

pretty.table.numbers = function(d)
   {d = copy(d)
    for (col in names(d)) d[, (col) :=
       (if (is.integer(get(col)))
            scales::comma(get(col), accuracy = 1)
        else if (str_detect(col, "\\ABias"))
            sprintf("%+.03f", get(col))
        else if (is.double(get(col)))
            sprintf("%.03f", get(col))
        else
            get(col))]
    d}

read_satellite_raster = function(
        satellite.product, tile, path, overpass = NA_integer_)
   {r1 = terra::rast(path)

    if (multipass.sat(satellite.product))
       {if (is.na(overpass))
          # Use a default overpass.
            overpass = 1L
        if (overpass == 1 && !any(str_detect(names(r1), "_1\\Z")))
          # This file has only one overpass, so `overpass` needs to be
          # NA for the following logic.
            overpass = NA_integer_}

    if (!is.na(overpass))
      # Reduce to the layers for this overpass.
       {regex = paste0("_", overpass, "\\Z")
        r1 = r1[[str_subset(names(r1), regex)]]
        names(r1) = str_remove(names(r1), regex)}

    if (satellite.product == "geonexl2")
      # Rename a layer for consistency with MCD19A2.
        names(r1) = str_replace(names(r1), "AOT_Uncertainty", "AOD_Uncertainty")

    # Get the other (lower-resolution) layers.
    r2 = terra::rast(lapply(
        c("cosSZA", "cosVZA", "RelAZ", "Scattering_Angle", "Glint_Angle"),
        function(vname)
           {r = terra::rast(paste0(
                str_replace(r1@ptr$filenames()[1], ":grid1km:[A-Za-z0-9_]+\\Z",
                    ":grid5km:"),
                vname))
            if (!is.na(overpass))
                r = r[[overpass]]
            `names<-`(r, vname)}))

    # Combine all the layers.
    fix.coordinates = function(r)
       {if (satellite.product != "geonexl2")
            return(r)
        # The coordinates are all wrong. Fix them.
        # See "Example for 0.01x0.01 (1km) grid" at
        # https://web.archive.org/web/20220119040921/https://www.nasa.gov/geonex/dataproducts
        terra::crs(r) = paste0("epsg:", crs.lonlat)
        h0 = -180
        v0 = 60
        size = 6
        m = stringr::str_match(tile, "h(\\d+)v(\\d+)")
        assert(!anyNA(m))
        h = as.integer(m[,2])
        v = as.integer(m[,3])
        terra::ext(r) = c(
            h0 + (h + c(0, 1))*size,
            v0 - (v + c(1, 0))*size)
        r}
    c(
       fix.coordinates(r1),
       terra::disagg(fix.coordinates(r2), nrow(r1) / nrow(r2)))}

#' Compute the grid on which all predictions will be made.
make_pred_grid = function(satellite.product, earthdata.rows)
   {r = do.call(terra::merge, lapply(1 : nrow(earthdata.rows),
        function(i) with(earthdata.rows[i],
           {r = read_satellite_raster(satellite.product, tile, path)[[1]]
            r$tile = as.integer(tile)
            r$cell.local = seq_len(terra::ncell(r))
            r[[c("tile", "cell.local")]]})))
    r$tile = factor(as.integer(drop(r$tile[])),
        labels = levels(earthdata.rows$tile))
    r}

#' Train a full model using DART and 2-fold CV to select hyperparameters from
#' maximin Latin hypercube sample
#'
dart_full <- function(
  data_train,
  y_var,
  features,
  n_rounds = 100,
  progress = TRUE
){
  xgb_threads <- get.threads()
  data_train = data_train[, c(features, y_var), with = F]
  simplify.dt.for.xgboost(data_train)
  xdc_out <- xgboost.dart.cvtune(
    # by default, gives 100 rounds, and it is enough by experience
    n.rounds = n_rounds,
    d = data_train,
    dv = y_var,
    ivs = features,
    progress = progress,
    nthread = xgb_threads)
  # record the prediction
  preds = xdc_out$pred.fun(data_train)
  #rmse_full <- sqrt(mean((preds - data_train[[y_var]])^2))

  # SHAP
  shap_pred <- as.data.table(xdc_out$pred.fun(data_train, predcontrib = TRUE, approxcontrib = FALSE))
  shap_bias <- first(shap_pred$BIAS)
  shap_pred[, BIAS := NULL]

  # features ranked by SHAP
  mean_shaps = colMeans(abs(shap_pred))
  ## feature names only:
  # features_rank_full_model <- names(mean_shaps)[order(mean_shaps, decreasing = TRUE)]
  # named vector:
  features_rank_full_model = mean_shaps[order(mean_shaps, decreasing = T)]

  list(features_rank_full_model = features_rank_full_model,
       y_preds = preds,
       shap_pred = shap_pred,
       shap_bias = shap_bias,
       model = xgb.save.raw(raw_format = "ubj", xdc_out$model))}

get_aoi_buffer = function(region)
   {buffer.size.m = 5000

    region = (
        if (region == "conus")
            get_conus()
        else if (startsWith(region, "/"))
            read_sf(region)
        else
           {path = tempfile(fileext = paste0(".",
                tools::file_ext(region)))
            assert(0 == download.file(region, path))
            read_sf(path)})

    # Unify the region and buffer it out.
    st_buffer(dist = buffer.size.m, st_union(region))}

#' get_conus
#' @return sf object of states in CONUS
get_conus <- function(){

  x = read_sf(paste0("/vsizip/", download(
      "https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_us_state_20m.zip",
        # Linked to from https://www.census.gov/geographies/mapping-files/time-series/geo/cartographic-boundary.html
      "conus.zip")))
      x[!(x$STUSPS %in% c("AK", "HI", "PR")),]
}

in.sf = function(x, y, crs, region)
# Return a logical vector indicating whether each (x, y) point is in the
# specified region (an `sf` object).
    drop(st_intersects(sparse = F,
        convert.crs(cbind(x, y), crs, terra::crs(region), sf = T),
        region))

simplify.dt.for.xgboost = function(d)
    if ("time.sat" %in% colnames(d))
        d[, time.sat := as.double(time.sat)]
