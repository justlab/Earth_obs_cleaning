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
  stations = stations[Site_Name != "CART_SITE"]
    # This station is in the same place as "Cart_Site" and has
    # relatively few observations. It's not clear to what degree it's
    # a duplicate.
  aerpts = st_transform(crs = st_crs(reg_polygon),
      st_as_sf(stations, coords = c("lon", "lat"), crs = crs.lonlat))
  aerpts <- st_transform(crs = new_crs,
      aerpts[st_intersects(aerpts, reg_polygon, sparse = FALSE),])
  `rownames<-`(aerpts[withr::with_collate("C", order(aerpts$Site_Name)),],
      NULL)
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

#' Return a matrix classifying cells as being within a distance from the center.
#'
#' Assumes square matrix with odd number of rows & columns. The width (and
#' height) will be 2x the radius plus one cell.
#'
#' @param radius distance, in kilometers, from center cell to classify as inside
#'   circle
#' @param cellsize dimensions of raster cells in meters. Pixels assumed square.
#' @param matrix if TRUE, return a logical matrix. If FALSE (default), return a
#'   data.table with two columns.
#' @return a data table with offsets from the central cell that are within the
#'   specified radius
circle_mat = function(radius, cellsize = 926.6254, matrix = FALSE){
  radius = radius * 1000
  width = ceiling(radius/cellsize*2+1)
  if(width%%2 == 0) width <- width + 1 # ensure odd row and column count
  height = width # assuming square matrix
  circle_df = expand.grid(mrow = 1:height, mcol = 1:width)
  center_val = median(circle_df$mrow) # assuming square matrix
  circle_df$dist_cells = sqrt((circle_df$mrow - center_val)^2 + (circle_df$mcol - center_val)^2)
  circle_df$circ = ifelse(circle_df$dist_cells * cellsize <= radius, TRUE, FALSE)

  if(matrix == TRUE){
    circle_mat = matrix(as.numeric(circle_df$circ), height, width)
    circle_mat
  } else {
    setDT(circle_df)
    circle_df[, offset_x := mcol - center_val]
    circle_df[, offset_y := mrow - center_val]
    circle_df[circ == TRUE, .(offset_x, offset_y)]
    circle_df
  }
}

# Extent to create a raster version of focal filter matrix
# May extend beyond the AOD raster
get_focal_extent = function(x, c1, r1, radw){
  center_x = xFromCol(x, c1)
  center_y = yFromRow(x, r1)
  r = res(x)
  xn <- center_x - ((0.5 * r[1]) + (r[1] * radw))
  xx <- center_x + ((0.5 * r[1]) + (r[1] * radw))
  yn <- center_y - ((0.5 * r[2]) + (r[2] * radw))
  yx <- center_y + ((0.5 * r[2]) + (r[2] * radw))
  ext(c(sort(c(xn, xx))), sort(c(yn, yx)))
}

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
        .(ground$site), on = "Site_Name",
        as.data.table(st_coordinates(geometry))]]
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

    message("Reading predictors from satellite data")
    d = rbindlist(pblapply(
        cl = n.workers,
        split(d, by = c("sat.files.ix", "overpass")),
        function(chunk) chunk
            [, c("y.sat", vnames, "y.sat.mean", "y.sat.present") :=
               {r = read_satellite_raster(
                    satellite.product,
                    satellite_hdf_files[chunk$sat.files.ix[1], tile],
                    satellite_hdf_files[chunk$sat.files.ix[1], path],
                    overpass[1])
                cbind(
                    r[[c(y.sat.name, vnames)]][cell.local],
                    rbindlist(lapply(cell.local, function(cell)
                       {rc = terra::rowColFromCell(r, cell)
                        row = rc[,1] + (-window.radius : window.radius)
                        col = rc[,2] + (-window.radius : window.radius)
                        values = r[[y.sat.name]][
                            row[0 <= row & row <= nrow(r)],
                            col[0 <= col & col <= ncol(r)]][[1]]
                        data.frame(
                            y.sat.mean = mean(values, na.rm = T),
                            y.sat.present = mean(!is.na(values)))})))}]
            [!is.na(y.sat)]))
              # Keep only cases where we have the satellite
              # outcome of interest.

    if ("time.diff" %in% colnames(d))
        d[, time.diff := NULL]

    if ("AOD_QA" %in% colnames(d))
       {d[, qa_best := bitwAnd(AOD_QA,
            bitwShiftL(strtoi("1111", base = 2), 8)) == 0]
              # Page 13 of https://web.archive.org/web/20200927141823/https://lpdaac.usgs.gov/documents/110/MCD19_User_Guide_V6.pdf
        d[, AOD_QA := NULL]}

    d}

derive_mcd19_vars = function(aer_data, n.workers, ...)
  {aer_data[, chunk := match(aer_date, sort(unique(aer_date))) %% n.workers]
   d = rbindlist(Filter(nrow, parallel::mclapply(mc.cores = n.workers,
       split(aer_data, by = "chunk"),
       function(chunk)
           chunk[, by = aer_date, .SDcols = colnames(chunk),
               derive_mcd19_vars_1day(.SD, ...)])))
   assert(nrow(unique(d[, .(Site_Name, stn_time)])) == nrow(d))
   assert(nrow(unique(d[, .(cell, overpass_time)])) == nrow(d))
   d}

#' Roll join a single day of AERONET data to nearest MCD19A2 overpass and
#' calculate derived values from MCD19A2 within specified distances from the
#' AERONET station.
#'
#' @param aer_data A single day of AERONET observations
#' @param load_sat input "terra" or "aqua"
#' @param buffers_km vector of buffer radii in kilometers
#' @param aer_stn data table of AERONET station names and reference unique ID
#'   for satellite AOD cells
#' @param satellite_hdf_files data table returned by `get_earthdata`
#' @param agg_level how much to aggregate the input MODIS AOD to use for large
#'   focal radii
#' @param agg_thresh use the aggregated AOD when `(radius * 1000 / agg_level / mcd_res) > agg_thresh`,
#'   where `radius` is in km, and `mcd_res` (resolution of input AOD) is
#'   in meters.
#' @param vrt_path directory to store overpass VRTs that reference HDF files.
#' @param rolldiff_limit maximum time between the AERONET and satellite
#'   observations
derive_mcd19_vars_1day = function(aer_data, load_sat, buffers_km, aer_stn, satellite_hdf_files,
                                  agg_level, agg_thresh, vrt_path,
                                  rolldiff_limit = as.difftime(7.5, units = 'mins')){
  set.seed.with.obj(list(
      "derive_mcd19_vars_1day", aer_data$aer_date[1], load_sat))
  # 1. Prepare AERONET data
  aer_stn = st_sf(aer_stn) # testing whether passing in a DT version avoids error with vctrs package
  sv_aer = vect(aer_stn)
  setDT(aer_data)
  if(nrow(aer_data) == 0){
    return(data.table())
  }
  aer_data = interpolate_aod(aer_data, aer_stn)

  # Get the date of this chunk of AERONET data, find HDFs from same date
  if(aer_data[, uniqueN(aer_date)] > 1) stop('Use tar_group_by and pattern=map to send one date at a time to derive_mcd19_vars')
  this_date = aer_data[1, aer_date]

  hdf_files = (satellite_hdf_files
      [.("terra.and.aqua", this_date)]
      [file.size(path) != 0])

  if(nrow(hdf_files) == 0){
    return(data.table())
  }
  # 2. Join AERONET to satellite AOD
  binDT = bin_overpasses(hdf_files)
  binDT = binDT[sat == load_sat]
  if(nrow(binDT) == 0){
    return(data.table())
  }

  day_op = get_overpasses_vrts(hdf_files$path, binDT, load_sat, vrt_path)

  # lapply by overpass, which are the groups in day_op
  tryCatch({
  day_rasters = mapply(FUN = rasterize_vrts,
                       op_vrts = day_op,
                       opDT = lapply(1:length(day_op), function(i) binDT[overpass_bin == i]),
                       SIMPLIFY = FALSE)
  }, error = function(e) {
    emsg = paste0('this_date = ', this_date, '\n',
           'length(day_op) = ', length(day_op), '\n',
           'binDT[, uniqueN(overpass_bin)] = ', binDT[, uniqueN(overpass_bin)], '\n',
           'error = ', e)
    stop(emsg)
  })
  rm(day_op)

  # roll join will update value of the time column in X to the value of the time
  # column in i. Copy AERONET's stn_time to a column with the same name as
  # MCD19's time column so we can compare the time difference afterwards using
  # meaningful column names.
  setnames(aer_data, 'stn_time', 'overpass_time')
  aer_data[, stn_time := overpass_time]

  # extract AOD overpass values by each overpass
  overpass_stats <- function(single_op){

    # extract values from each raster for the overpass
    op_vals = lapply(single_op, function(op){
      # only record cell ID for the raster stack with AOD
      if(length(grep('Optical_Depth', names(op))) > 0){
        dt = setDT(terra::extract(op, sv_aer, cells = TRUE))
      } else {
        dt = setDT(terra::extract(op, sv_aer, cells = FALSE))
      }
      setkey(dt, ID)
    })
    op_vals = Reduce(merge.data.table, op_vals)
    op_vals[, Site_Name := aer_stn$Site_Name]
    setnames(op_vals, 'Optical_Depth_047', 'MCD19_AOD_470nm')
    op_vals = op_vals[!is.na(MCD19_AOD_470nm)]
    op_vals[, overpass_time := as.POSIXct(overpass_time, tz = 'UTC', origin = '1970-01-01')]

    # roll join AERONET to satellite AOD
    setkey(op_vals, Site_Name, overpass_time)
    rj <- aer_data[op_vals, roll = 'nearest', nomatch = 0] # using keys to join

    # Only keep roll joins within specified time difference limit
    rj[, rj_difftime := overpass_time - stn_time]
    rj = rj[abs(rj_difftime) <= rolldiff_limit, ]

    # Only keep one case per satellite observation. Duplicates occur
    # when there are multiple stations in one cell of the satellite
    # grid.
    rj = rj[order(abs(rj_difftime), Site_Name),
        by = .(overpass_time, cell), head(.SD, 1)]

    if(nrow(rj) > 0){
      # 3. Derived values from focal statistics on satellite AOD

      # find the AOD layer in the list of rasters for the overpass
      mcd_aod = single_op[unlist(lapply(single_op,
        function(s) 'Optical_Depth_047' %in% names(s)))][[1]][['Optical_Depth_047']]
      mcd_nna = !is.na(mcd_aod)
      mcd_res = res(mcd_aod)[1] # assuming square pixels
      # Aggregate AOD raster for use in larger focal radii
      mcd_agg_aod = terra::aggregate(mcd_aod, fact = agg_level, fun = 'mean', na.rm = TRUE)
      mcd_agg_nna = terra::aggregate(mcd_nna, fact = agg_level, fun = 'mean', na.rm = TRUE)

      filter_matrices = list()
      use_agg = list()
      for(radius in buffers_km){
        if(radius * 1000 / agg_level / mcd_res < agg_thresh){
          # use original resolution
          use_agg[[paste0('r', radius)]] <- FALSE
          this_mat = circle_mat(radius, cellsize = mcd_res, matrix = TRUE)
        } else {
          # use aggregated resolution
          use_agg[[paste0('r', radius)]] <- TRUE
          this_mat = circle_mat(radius, cellsize = mcd_res * agg_level, matrix = TRUE)
        }
        this_mat[this_mat == 0] <- NA
        filter_matrices[[paste0('r', radius)]] <- this_mat
      }
      rm(radius)

      # Calculate focal statistics for each unique cell containing an AERONET station
      # returns a data.table to merge (send unique values in case any stations in same cell)
      station_focal <- function(cellid){
        agg_cellid = cellFromXY(mcd_agg_aod, xyFromCell(mcd_aod, cellid))
        focal_list = list(cellid = cellid, agg_cellid = agg_cellid)
        # buffers should be in km, resolution should be in m
        for(radius in buffers_km){
          filter_mat = filter_matrices[[paste0('r', radius)]]
          focw = ncol(filter_mat)  # width of focal matrix (assumed square matrix)
          radw = (focw - 1) / 2    # width of the radius in cells, without central cell
          # Extract matrix of AOD values around station
          if(use_agg[[paste0('r', radius)]] == FALSE){
            # original AOD raster resolution
            crid = rowFromCell(mcd_aod, cellid) # central row id
            ccid = colFromCell(mcd_aod, cellid) # central column id
            filter_raster = rast(filter_mat, crs = crs_sinu,
                                 extent = get_focal_extent(mcd_aod, ccid, crid, radw))
            aod_extract = crop(mcd_aod, filter_raster)
            nna_extract = crop(mcd_nna, filter_raster)
          } else {
            # aggregated AOD raster
            crid = rowFromCell(mcd_agg_aod, agg_cellid) # central row id
            ccid = colFromCell(mcd_agg_aod, agg_cellid) # central column id
            filter_raster = rast(filter_mat, crs = crs_sinu,
                                 extent = get_focal_extent(mcd_agg_aod, ccid, crid, radw))
            aod_extract = crop(mcd_agg_aod, filter_raster)
            nna_extract = crop(mcd_agg_nna, filter_raster)
          }
          filter_crop = crop(filter_raster, aod_extract)
          mult_aod = aod_extract * filter_crop
          mult_nna = nna_extract * filter_crop
          focal_list[[paste0('Mean_AOD', radius, 'km')]] <- mean(mult_aod[], na.rm = TRUE)
          focal_list[[paste0('pNonNAAOD', radius, 'km')]] <- mean(mult_nna[], na.rm = TRUE)
        }
        setDT(focal_list)
        focal_list
      }
      focal_stats = rbindlist(lapply(rj[, unique(cell)], station_focal))
      setkey(focal_stats, cellid)
      setkey(rj, cell)
      rj = rj[focal_stats]

      # difference between central cell satellite AOD and mean AOD in buffers
      for(distx in buffers_km){
        rj[, c(paste0('diff_AOD', distx, 'km')) := MCD19_AOD_470nm - get(paste0('Mean_AOD', distx, 'km'))]
      }
      rj[, c('ID', 'x_satcrs', 'y_satcrs') := NULL] # remove AERONET station row index and coordinates
      rj
    } else {
      # This will create an extra "V1" column in rowbound final target with all NAs,
      # but it gets around the error of writing NULL to FST format
      data.table()
    }
  }
  rbindlist(lapply(day_rasters, overpass_stats), fill = TRUE)
}

create_qc_vars <- function(dt){
  dt[, qa_bits := bitwAnd(AOD_QA, strtoi("111100000000", base = 2))]
  dt[, qa_best := 0]
  dt[qa_bits == 0 , qa_best := 1]

  dt[, qa_bits := bitwAnd(AOD_QA, strtoi("11000", base = 2))]
  dt[qa_bits==0, qa_lwsi := "land"]
  dt[qa_bits==8, qa_lwsi := "water"]
  dt[qa_bits==16, qa_lwsi := "snow"]
  dt[qa_bits==24, qa_lwsi := "ice"]
  return(dt)
}

#' Calculate DV, remove rows without DV, and calculate dayint and QC var IVs.
#'
#' Optionally subset the MCD19A2 observations to those with dates in the
#' \code{date_range} vector.
#'
#' @param dt data.table of MCD19A2 observations
#' @param date_range vector of Date type containing all dates to retain in the
#'   subset
#' @return data.table with DV and additional IV columns ready for running CV
prepare_dt <- function(dt, date_range = NULL){
  setnames(dt, "Optical_Depth_047", "MCD19_AOD_470nm", skip_absent=TRUE)
  if(!is.null(date_range)){
    dt = dt[aer_date %in% date_range, ]
  }
  # The dependent variable: diff_AOD = MCD19 - AERONET = Optical_Depth_047 - AOD_470nm
  dt[, diff_AOD := MCD19_AOD_470nm - AOD_470nm]
  dt <- create_qc_vars(dt)
  dt[, dayint := as.integer(as.Date(overpass_time))]
  dt = dt[!is.na(get(y_var))]
  # drop overpass matchups for pixels classified as water
  dt = dt[qa_lwsi != "water"]
  dt
}

#' Do initial CV with hyperparameter selection with DART and rank features.
#'
#' Must provide either stn_var or day_var for binning folds.
#'
#' @param data data.frame of observations with DV and IVs
#' @param y_var column name to predict
#' @param features character vector of column names to train on
#' @param stn_var if binning by stations, provide the station column name
#' @param day_var if binning by days, provide the day column name
#' @param absolute whether to use absolute loss (rather than square loss)
# Depends on functions in xgboost_cv_RFE.R
initial_cv_dart <- function(
  data,
  y_var,
  features,
  k_fold = 5,
  n_rounds = 100,
  stn_var = NULL,
  day_var = NULL,
  absolute = FALSE,
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
                              absolute = absolute,
                              seed = (if (absolute) 400 else 1234),
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

# Prediction ####

#' #' Given an input SF object, return a terra:::SpatExtent in a new CRS,
#' #' optionally extended in all directions by `expand_distance`, using the units
#' #' of the `new_crs`.
#' get_aoi_ext <- function(shape, new_crs, expand_distance = 0){
#'   shape = st_transform(shape, new_crs)
#'   bbox = st_bbox(shape)
#'   if(expand_distance != 0){
#'     mod_box = st_bbox(c(xmin = -expand_distance, xmax = expand_distance,
#'                         ymin = -expand_distance, ymax = expand_distance))
#'     bbox = bbox + mod_box
#'   }
#'   ext(bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax)
#' }

#' Given an input SF object, return a terra:::SpatExtent. Optionally reproject
#' to `to_crs`.
sf_to_ext <- function(sf, to_crs = NULL){
  if(!is.null(to_crs)) sf = st_transform(sf, st_crs(to_crs))
  sv = vect(sf)
  ext(sv)
}

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
                str_replace(r1@ptr$filenames[1], ":grid1km:[A-Za-z0-9_]+\\Z",
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
    # We need to call `as.integer` again: https://github.com/rspatial/terra/issues/866
    r$cell.local = as.integer(r$cell.local[])
    r$tile = factor(as.integer(drop(r$tile[])),
        labels = levels(earthdata.rows$tile))
    r}

#' Prepare prediction table. Will predict values for every overpass in the
#' specified region and dates.
#'
#' @param features character vector of column names to train on
#' @param buffers_km vector of buffer radii in kilometers
#' @param satellite_hdf_files data table returned by `get_earthdata`
#' @param vrt_path directory to store overpass VRTs that reference HDF files.
#' @param this_date a date to make predictions for
#' @param sat input "terra" or "aqua"
#' @param agg_level how much to aggregate the input MODIS AOD to use for large
#'   focal radii
#' @param agg_thresh use the aggregated AOD when `(radius * 1000 / agg_level /
#'   mcd_res) > agg_thresh`, where `radius` is in km, and `mcd_res` (resolution
#'   of input AOD) is in meters.
#' @param aoi sf shape of the region the model was trained over.
#' @param pred_bbox SpatExtent of the region to predict. If NULL, will crop
#'   using the shape provided in `aoi`.
pred_inputs <- function(features, buffers_km, satellite_hdf_files, vrt_path, load_sat,
                        this_date, agg_level, agg_thresh, aoi,
                        pred_bbox = NULL){

  # Get overpasses for the date
  if(length(this_date)>1) stop('Unexpected this_dates vector longer than 1')

  hdf_files = (satellite_hdf_files
      [.("terra.and.aqua", this_date)]
      [file.size(path) != 0])

  if(nrow(hdf_files) == 0){
    return(data.table())
  }

  binDT = bin_overpasses(hdf_files)
  binDT = binDT[sat == load_sat]
  if (!nrow(binDT))
      return(data.table())
  day_op = get_overpasses_vrts(hdf_files$path, binDT, load_sat, vrt_path)

  # lapply by overpass, which are the groups in day_op
  tryCatch({
    day_rasters = mapply(FUN = rasterize_vrts,
                         op_vrts = day_op,
                         opDT = lapply(1:length(day_op), function(i) binDT[overpass_bin == i]),
                         SIMPLIFY = FALSE)
  }, error = function(e) {
    emsg = paste0('this_date = ', this_date, '\n',
                  'length(day_op) = ', length(day_op), '\n',
                  'binDT[, uniqueN(overpass_bin)] = ', binDT[, uniqueN(overpass_bin)], '\n',
                  'error = ', e)
    stop(emsg)
  })
  rm(day_op)

  aoi = sf_to_ext(aoi, crs_sinu)
  if(is.null(pred_bbox)) pred_bbox = aoi

  terraOptions(progress = 0)

  # Prepare predictors for each overpass
  op_to_table <- function(op_id){
    raslist = day_rasters[[op_id]]
    # note: hardcoded list item names (based on MCD19A2 resolution) and harcoded
    #   disaggregation factors
    r_relaz = terra::disagg(raslist$r4633, 5)
    r_times = terra::disagg(raslist$time, 1200)
    rstack = rast(list(raslist$r926, r_relaz, r_times))
    rm(r_relaz, r_times)
    # expand bbox to accomodate largest focal filter width
    # SpatExtent automatically subtracts the value on the minimum sides of the box, unlike st_bbox
    exp_bbox = pred_bbox + max(buffers_km)*1000
    rstack = terra::crop(rstack, exp_bbox)
    # keep a reference to the AOD layer for use in focal stats
    ras_aod = rstack$Optical_Depth_047
    aod_res = res(ras_aod)[1] # assuming square pixels
    # non-NA AOD layer
    ras_nna = !is.na(ras_aod)
    # aggregate AOD and non-NA for focal use
    agg_aod = terra::aggregate(ras_aod, fact = agg_level, fun = 'mean', na.rm = TRUE)
    agg_nna = terra::aggregate(ras_nna, fact = agg_level, fun = 'mean', na.rm = TRUE)
    # focal statistics
    focal_list = list()
    for(radius in buffers_km){
      if(radius * 1000 / agg_level / aod_res < agg_thresh){
        # use original resolution
        this_mat = circle_mat(radius, cellsize = aod_res, matrix = TRUE)
        foc_aod <- terra::focal(ras_aod, this_mat, fun = 'mean', na.rm = TRUE)
        foc_nna <- terra::focal(ras_nna, this_mat, fun = 'mean', na.rm = TRUE)
        focal_list[[paste0('Mean_AOD', radius, 'km')]]  <- foc_aod
        focal_list[[paste0('pNonNAAOD', radius, 'km')]] <- foc_nna
      } else {
        # use aggregated resolution
        this_mat = circle_mat(radius, cellsize = aod_res * agg_level, matrix = TRUE)
        foc_aod_agg <- terra::focal(agg_aod, this_mat, fun = 'mean', na.rm = TRUE)
        foc_nna_agg <- terra::focal(agg_nna, this_mat, fun = 'mean', na.rm = TRUE)
        fog_aod = terra::disagg(foc_aod_agg, agg_level)
        fog_nna = terra::disagg(foc_nna_agg, agg_level)
        rm(foc_nna_agg, foc_aod_agg)
        focal_list[[paste0('Mean_AOD', radius, 'km')]]  <- foc_aod
        focal_list[[paste0('pNonNAAOD', radius, 'km')]] <- foc_nna
      }
    }
    rm(radius)
    # stack focal results with the expanded crop above
    rstack = rast(list(rstack, rast(focal_list)))

    # crop to actual AOI
    rstack = terra::crop(rstack, pred_bbox)

    # convert to data.table
    rasDT = setDT(as.data.frame(rstack, na.rm = FALSE, xy = TRUE))
    rasDT = rasDT[!is.na(Optical_Depth_047)]
    setnames(rasDT, 'Optical_Depth_047', 'MCD19_AOD_470nm')
    rasDT[, overpass_time := as.POSIXct(overpass_time, tz = 'UTC', origin = '1970-01-01')]
    rasDT[, dayint := as.integer(as.Date(overpass_time))]

    # calc diffAOD for each buffer
    for(radius in buffers_km){
      rasDT[, c(paste0('diff_AOD', radius, 'km')) := MCD19_AOD_470nm - get(paste0('Mean_AOD', radius, 'km'))]
    }

    # convert AOD_QA to qa_best
    create_qc_vars(rasDT)

    if (nrow(rasDT))
        # overpass bin
        rasDT[, op_id := op_id]
    rasDT
  }
  rbindlist(fill = TRUE, lapply(1:length(day_rasters), op_to_table))
}

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

# Maps ####

#' Compare adjusted AOD to original
#' @param data the data.table output of prediction
#' @param viz_op the single overpass to visualize from the selected date
#' @return vertically stacked ggplots comparing original and adjusted MCD19A2 AOD
ggplot_orig_vs_adj = function(data, viz_op, grid){
  data = data[overpass == viz_op,
              .(cell, MCD19_AOD_470nm = value_old * 0.00001, MCD19_adjust = value_new * 0.00001)]
  data[, c("x", "y") := data.table(terra::xyFromCell(grid, cell))]
  orig = simple.pred.map(data, fillvar = 'MCD19_AOD_470nm')
  adj = simple.pred.map(data, fillvar = 'MCD19_adjust')
  title = ggdraw() + draw_label(paste('Overpass', viz_op))
  cowplot::plot_grid(title, orig, adj, ncol = 1, align = 'v', axis = 'r')
}

# adapted from CONUS_air:plots.R
simple.pred.map = function(preds, fillvar, xvar = 'x', yvar = 'y',
                           scale.args = list()){
  ggplot(preds) +
  geom_raster(aes_string(xvar, yvar, fill = fillvar)) +
  do.call(scale_fill_distiller,
          c(list(palette = "Spectral"), scale.args)) +
  coord_equal() +
  theme_void() +
  theme(legend.justification = "left")
}

#' Interactively compare adjusted AOD to original
#'
#' @param data the data.table output of prediction
#' @param viz_op the single overpass to visualize from the selected date
#' @param grid the base grid used for looking up cell coordinates
#' @param maxpixels resample raster version of predictions to approxmimately
#'   this many pixels and force the display in mapshot of this resolution
#' @return mapview object with a layer for the original and adjusted MCD19A2 AOD
mapshot_orig_vs_adj = function(data, viz_op, grid,
                               use_jenks = FALSE, maxpixels = NULL){
  data = data[overpass == viz_op,
              .(cell, MCD19_AOD_470nm = value_old * 0.00001, MCD19_adjust = value_new * 0.00001)]
  data[, c("x", "y") := data.table(terra::xyFromCell(grid, cell))]
  data[, cell := NULL]
  setcolorder(data, c("x", "y"))
  ras = rasterFromXYZ(data, crs = crs_sinu)

  # resample raster
  if(!is.null(maxpixels)){
    mapviewOptions(georaster = TRUE, mapview.maxpixels = maxpixels)
    ras = sampleRegular(ras, maxpixels, asRaster = TRUE)
  }

  ub = unified_breaks(ras, n_classes = 10, viridis::inferno, use_jenks = use_jenks)
  get_mapview = function(i){
    mapview::mapview(ras[[i]], na.color = '#AAAAAA00', alpha.regions = 1,
                     col.regions = ub[[i]]$colors,
                     at = ub[[i]]$breaks, layer.name = names(ras[[i]]))}
  maps = lapply(1:2, get_mapview)
  if(!dir.exists(intermediate.path('mapshot'))) dir.create(intermediate.path('mapshot'), recursive = T)
  mapshot_path = intermediate.path('mapshot',
                      paste0('orig_adj_MCD19_',
                             targets:::digest_obj64(list(data, use_jenks, maxpixels)),
                             '.html'))
  mapshot(maps[[1]] + maps[[2]], url = mapshot_path)
  out_size_MB = file.size(mapshot_path)/1024^2
  if(out_size_MB>500) warning('The output file', mapshot_path, 'is', round(out_size_MB,0), 'MiB in size')
  mapshot_path
}

#' Get unified breaks and color scheme for layer ranges that partially overlap
#' @param ras RasterStack or RasterBrick to prepare a unified set of class
#'   breaks for all its layers
#' @param n_classes number of classes, a numeric of length 1
#' @param color_func function to use to calculate colors, should take a single
#'   number as input
#' @return list with length equal to count of raster layers. Each item contains
#'   \code{$breaks}: numeric vector of class breaks, and \code{$colors} a
#'   character vector of hex-encoded colors.
unified_breaks = function(ras, n_classes, color_func, use_jenks = FALSE){
  rrange = range(ras[], na.rm = TRUE)
  if(use_jenks == TRUE){
    if(n_classes < 2) stop('cannot use n_classes < 2 with use_jenks = TRUE')
    # get all raster data into a single vector, removing NA values
    datavec = base::Reduce(c, as.data.table(ras[])[!is.na(get(names(ras)[1]))])
    nb1 = rgeoda::natural_breaks(k = n_classes, df = as.data.table(datavec))
    full_breaks = c(min(datavec), nb1, max(datavec))
    rm(datavec)
  } else {
    full_breaks = seq(rrange[1], rrange[2],
                      (rrange[2]-rrange[1])/n_classes)
  }
  colors = color_func(n_classes)

  breaks_by_layer = function(layernum, ras, full_breaks, colors){
    layer_min = min(ras[[layernum]][], na.rm = TRUE)
    layer_max = max(ras[[layernum]][], na.rm = TRUE)
    lyrI = intersect(which(full_breaks >= layer_min), which(full_breaks <= layer_max))
    lyrB = full_breaks[lyrI]

    if(full_breaks[min(lyrI)] > layer_min){
      lyrB = c(layer_min, lyrB)
      lyrI = c(min(lyrI)-1, lyrI)
    }
    if(full_breaks[max(lyrI)] < layer_max){
      lyrB = c(lyrB, layer_max)
      lyrI = c(lyrI, max(lyrI) + 1)
    }
    list(breaks = lyrB, colors = colors[lyrI[1:(length(lyrI)-1)]])
  }

  lapply(1:nlayers(ras), FUN = breaks_by_layer,
         ras = ras, full_breaks = full_breaks, colors)
}

#' get_aoi_buffer
#' @param aoiname the region to to use for selection buffer
#' @return sf object of region buffer
get_aoi_buffer <- function(aoiname){
  switch(aoiname,
    "conus" = buff <- get_conus_buff(),
    stop("Unsupported aoi name:", aoiname)
  )
  buff
}

#' get_conus
#' @return sf object of states in CONUS
get_conus <- function(){

  x = read_sf(paste0("/vsizip/", download(
      "https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_us_state_20m.zip",
        # Linked to from https://www.census.gov/geographies/mapping-files/time-series/geo/cartographic-boundary.html
      "conus.zip")))
      x[!(x$STUSPS %in% c("AK", "HI", "PR")),]
}


#' get_conus_buff
#' @return sf object of buffered CONUS area
get_conus_buff <- function(){
  buffer.size.m = 5000

  x = get_conus()
  st_buffer(dist = buffer.size.m,
            st_transform(crs = crs.us.atlas,
                         st_union(x)))
}

simplify.dt.for.xgboost = function(d)
    if ("time.sat" %in% colnames(d))
        d[, time.sat := as.double(time.sat)]
