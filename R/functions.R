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

#' prepare AERONET observation data for interpolation of AOD470nm
#'
#' interpolating to wavelength 470 nm, input as 0.47
#'
#' @param aer_data the dataset contains all the AOD measurements and exact wavelengths
#' @param aer_stns simple feature collection of station locations
#' @return dataset of aer_data with pred, also with `nearest_refgrid` next to Site_Name
interpolate_aod <- function(aer_data, aer_stns){
  keepcols = grep('1020|1640', names(aer_data), invert = TRUE)
  aer_data <- aer_data[, ..keepcols]
  aer_data[, obs:= .I]
  setkey(aer_data, obs)
  vars_wv <- grep("Exact_Wavelengths_of_AOD", names(aer_data), value = TRUE)
  vars_aod_sub <- grep("AOD_", names(aer_data), value = TRUE)
  aer_data$N_NAAOD <- rowSums(!is.na(aer_data[, ..vars_aod_sub]))
  # create long-format data
  d <- melt(aer_data[, c(vars_aod_sub, vars_wv, "obs"), with = FALSE],
            measure = list(vars_aod_sub, vars_wv), value.name = c("aod", "exact_wv"))
  # choose non-NA and >0 AOD
  d <- d[!is.na(aod) & aod > 0]
  # ids for wavelengths observed on one side of 470nm
  obs_os <- d[, by = obs, .(oneside = (max(exact_wv)-0.47)*(min(exact_wv)-0.47)>0)][oneside == TRUE, obs]
  # the predict part:
  d2 <- d[, by = obs,
     {m = lm(log(aod) ~ log(exact_wv) + I((log(exact_wv))^2))
      if (m$rank == length(coef(m)))
        # Only use the model if it's full-rank.
          list(AOD_470nm = exp(predict(m,
              data.frame(exact_wv = 0.47))))}] # notice to put 0.47 not 470
  setkey(d2, obs)
  aer_data_wPred <- d2[aer_data]
  # set NA: At least 4 observations wanted for exterpolation
  aer_data_wPred[N_NAAOD < 4 & obs %in% obs_os, AOD_470nm := NA]
  # remove unnecessary variables:
  aer_data_wPred <- aer_data_wPred[, -c(vars_aod_sub, vars_wv), with = FALSE]

  # join station locations to observations
  setnames(aer_data_wPred, "AERONET_Site_Name", "Site_Name")
  setkey(aer_data_wPred, Site_Name)
  aer_sites = as.data.table(aer_stns)
  aer_sites[, c('x_sinu', 'y_sinu') := as.data.table(st_coordinates(geometry))]
  aer_sites = aer_sites[, .(Site_Name, x_sinu, y_sinu)]
  setkey(aer_sites, Site_Name)
  aer_data_wPred <- aer_sites[aer_data_wPred]
  setkey(aer_data_wPred, Site_Name, stn_time)
  aer_data_wPred
}

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

derive_mcd19_vars = function(aer_data, n.workers, ...)
  {aer_data[, chunk := match(aer_date, sort(unique(aer_date))) %% n.workers]
   d = rbindlist(parallel::mclapply(mc.cores = n.workers,
       split(aer_data, by = "chunk"),
       function(chunk)
           chunk[, by = aer_date, .SDcols = colnames(chunk),
               derive_mcd19_vars_1day(.SD, ...)]))
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
      rj[, c('ID', 'x_sinu', 'y_sinu') := NULL] # remove AERONET station row index and coordinates
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

  # output
  c(list(by_var = by_var, bin_list = bin_list, mDT_wPred = mDT), cv_results)
}

#' Summarize CV statistics on all initial_cv objects
#'
#' @param cv_list list of all CV output objects; may include every year, both
#'   satellites, and both objective functions
#' @return data.table summarizing CV statistics for each loss and satellite
cv_summary <- function(cv_list){
  stats_list = vector(mode = "list", length = length(cv_list))
  difftimes_list = vector(mode = "list", length = length(cv_list))
  for(i in 1:length(cv_list)){
    cv = cv_list[[i]]
    dt = data.table::copy(cv$mDT_wPred)
    dt[, aodhat := MCD19_AOD_470nm - diff_AOD_pred]
    stats = performance_metrics(dt,
                                ground_truth = "AOD_470nm",
                                eo_raw = "MCD19_AOD_470nm",
                                eo_pred = "aodhat")
    stats$sat <- str_extract(names(cv_list)[[i]], 'terra|aqua')
    stats$loss <- str_match(names(cv_list)[[i]], 'cv_(l[12])_')[,2]
    stats_list[[i]] <- stats
    difftimes_list[[i]] <- round(summary(as.numeric(abs(cv$mDT_wPred$rj_difftime))),0)
  }
  # prediction summary stats
  statsDT = rbindlist(lapply(stats_list, as.data.table))
  # setcolorder(statsDT, c('sat', 'loss',
  #                      'MAE_uncorr', 'MAE_corr', 'MAE_change',
  #                      'rmse', 'MAD_mcd19', 'MAD_aodhat', 'MAD_change'))
  # difftime distribution
  difftimeDT = rbindlist(lapply(difftimes_list, function(x) as.list(x)))
  difftimeDT[, c('sat', 'loss') := statsDT[, .(sat, loss)]]
  setcolorder(difftimeDT, c('sat', 'loss'))

  setkey(statsDT, sat, loss)
  setkey(difftimeDT, sat, loss)
  list(stats = statsDT, difftimes = difftimeDT)
}

#' Calculate performance metrics on a data.table
#'
#' Comparative cross-validation performance metrics with raw and corrected values.
#' Takes a \code{data.table} and can be run across strata with a \code{by} statement.
performance_metrics <- function(dt, ground_truth, eo_raw, eo_pred, digits = 3){
  truth = quote(get(ground_truth))
  raw = quote(get(eo_raw))
  pred = quote(get(eo_pred))
  mae = function(v1, v2) mean(abs(v1 - v2))
  mad = function(v1) mean(abs(v1 - median(v1)))
  mse = function(v1, v2) mean((v1 - v2)^2)
  # rmse = function(v1, v2) sqrt(mean((v1 - v2)^2))
  mean_bias = function(v1, v2) mean(v1 - v2)
  r = function(x) round(x, digits)
  list(
    # MAE_raw = r(mae(dt[, eval(raw)], dt[, eval(truth)])),
    # MAE_pred   = r(mae(dt[, eval(pred)], dt[, eval(truth)])),
    # MAD_truth= r(mad(dt[, eval(truth)])),
    # MAD_pred = r(mad(dt[, eval(pred)])),

    RMSE_raw = r(sqrt(mse(v1 = dt[, eval(raw)], v2 = dt[, eval(truth)]))),
    RMSE_pred= r(sqrt(mse(v1 = dt[, eval(pred)], v2 = dt[, eval(truth)]))),
    SD_truth   = r(sd(dt[, eval(truth)])),
    SD_pred  = r(sd(dt[, eval(pred)])),
    SD_raw   = r(sd(dt[, eval(raw)])),
    # pct_of_raw_mse = round(100 * (mse(dt[, eval(pred)], dt[, eval(truth)]) / mse(dt[, eval(raw)], dt[, eval(truth)])), 1),
    bias_raw = r(mean_bias(dt[, eval(raw)], dt[, eval(truth)])),
    bias_pred = r(mean_bias(dt[, eval(pred)], dt[, eval(truth)])),
    r_raw = r(cor(dt[, eval(raw)], dt[, eval(truth)])),
    r_pred = r(cor(dt[, eval(pred)], dt[, eval(truth)])),
    stn_count  = uniqueN(dt$stn),
    train_N = dt[, .N],
    median_truth = r(median(dt[, eval(truth)]))
  )
}

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

#' Compute the grid on which all predictions will be made. We assume
#' that 1-m precision is sufficient and round the coordinates (both
#' here and when matching other rasters to this grid).
make_pred_grid = function(raster.paths)
   {d = rbindlist(lapply(raster.paths, function(p)
       {r = terra::rast(p)
        lapply(as.data.table(round(xyFromCell(r, 1 : ncell(r)))),
            as.integer)}))
    setkey(d, x, y)
    d}

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
  first_date = min(data_train$aer_date)
  data_train = data_train[, c(features, y_var), with = F]

  xdc_out <- xgboost.dart.cvtune(
    # by default, gives 100 rounds, and it is enough by experience
    n.rounds = n_rounds,
    d = data_train,
    dv = y_var,
    ivs = features,
    progress = progress,
    nthread = xgb_threads)

  # manually select and store some params
  param_dart <- c(xdc_out$model$params[c(1,2,4, 6:11)], nrounds = xdc_out$model$niter)

  # save the model
  # hash important input and output to create a unique name
  model_out_path = intermediate.path(paste0('full_model_dart_',
                            targets:::digest_obj64(list(data_train, features, y_var, param_dart)),
                            '.xgb'))
  xgboost::xgb.save(xdc_out$model, model_out_path)

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
       model_out_path = model_out_path,
       first_date = first_date)
}

#' Predict the difference between MCD19 and AERONET
run_preds = function(full_model, features, grid, round_digits, data){
  if (!nrow(data))
      return()

  data = data[!is.na(MCD19_AOD_470nm)]
  data[, pred_date := as.Date(dayint, '1970-01-01')]

  predvec = predict(
      xgboost::xgb.load(full_model$model_out_path),
      as.matrix(data[, features, with = FALSE]))
  # join predictions
  r = function(v) as.integer(round(10^round_digits * v))
  data = data[, .(
      pred_date,
      overpass = op_id,
      cell = grid[.(round(data$x), round(data$y)), which = T],
      value_old = r(MCD19_AOD_470nm),
      value_new = r(MCD19_AOD_470nm - predvec))]
  assert(!anyNA(data))
  setkey(data, pred_date, overpass, cell)
  data
}

# Maps ####

#' Compare adjusted AOD to original
#' @param data the data.table output of prediction
#' @param viz_date the single date to visualize from the output predictions
#' @param op_id the single overpass to visualize from the selected date
#' @return vertically stacked ggplots comparing original and adjusted MCD19A2 AOD
ggplot_orig_vs_adj = function(data, viz_date, viz_op){

  data = data[pred_date == viz_date & op_id == viz_op,
              .(x, y, MCD19_AOD_470nm, MCD19_adjust)]
  orig = simple.pred.map(data, fillvar = 'MCD19_AOD_470nm')
  adj = simple.pred.map(data, fillvar = 'MCD19_adjust')
  title = ggdraw() + draw_label(paste(as.character(viz_date), 'Overpass', viz_op))
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
#' @param viz_date the single date to visualize from the output predictions
#' @param op_id the single overpass to visualize from the selected date
#' @param maxpixels resample raster version of predictions to approxmimately
#'   this many pixels and force the display in mapshot of this resolution
#' @return mapview object with a layer for the original and adjusted MCD19A2 AOD
mapshot_orig_vs_adj = function(data, viz_date, viz_op,
                               use_jenks = FALSE, maxpixels = NULL){

  data = data[pred_date == viz_date & op_id == viz_op,
              .(x, y, MCD19_AOD_470nm, MCD19_adjust)]
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
