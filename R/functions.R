# Contains most functions required for the AOD cleaning targets workflow

#' Return AERONET points that intersect the polygon describing the area of
#' interest, and include the cell ID (\code{idM21pair0}) from the reference
#' CONUS MODIS grid.
#'
#' @param stations data.table including coordiantes (\code{lon, lat}) points of
#'   AERONET stations
#' @param reg_polygon SF polygon defining the area of interest
#' @param refgrid_path FST with coordinates of all raster cells in the area of
#'   interest
#' @param refras_path path to MODIS reference raster stack
#' @return a subset of geometry points within the given polygon
select_stations <- function(stations, reg_polygon, refgrid_path, refras_path){
  aerpts = st_as_sf(stations, coords = c("lon", "lat"), crs = 4326)
  if(st_crs(aerpts) != st_crs(reg_polygon)){ # polygons are in nad83
    aerpts <- st_transform(aerpts, crs = st_crs(reg_polygon))
  }
  aerpts <- aerpts[st_intersects(aerpts, reg_polygon, sparse = FALSE),]
  aerpts <- station_cell_ids(aerpts, refgrid_path, refras_path) # now in sinusoidal CRS
  aerpts
}

#' Get the unique ID of the grid cell an AERONET station is in.
#'
#' @param stations_sf SF points of AERONET stations
#' @param refgrid_path FST with coordinates of all raster cells in the area of
#'   interest
#' @param refras_path path to MODIS reference raster stack
#' @return SF points of AERONET stations with unique MODIS grid cell IDs added
#'   (\code{idM21pair0})
station_cell_ids = function(stations_sf, refgrid_path, refras_path){
  gras = raster(refras_path, band = 4)
  stations_sf = st_transform(stations_sf, crs = st_crs(gras))
  stations_sf$cell_index = raster::extract(gras, as(stations_sf, 'Spatial'))

  gDT = read_fst(refgrid_path, as.data.table = TRUE, columns = c('idM21pair0', 'cell_index'))
  setkey(gDT, cell_index)
  setDT(stations_sf)
  setkey(stations_sf, cell_index)
  stations_sf[gDT, idM21pair0 := idM21pair0]
  stations_sf[, cell_index := NULL]
  setkey(stations_sf, Site_Name)
  stations_sf = st_sf(stations_sf)
  # previous row names no longer meaningful after resorting
  attributes(stations_sf)$row.names <- 1:nrow(stations_sf)
  stations_sf
}

# most basic variables to read
vars0 <- c("Date(dd:mm:yyyy)", "Time(hh:mm:ss)", "Day_of_Year","AERONET_Site_Name",
           "Site_Elevation(m)", "Ozone(Dobson)", "NO2(Dobson)",
           "Solar_Zenith_Angle(Degrees)", "Precipitable_Water(cm)")
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
  stn_names = unique(stations$Site_Name)
  # open files containing names of stations in the specified region
  aer_files_dir <- sapply(paste0(unique(stn_names),".*\\.lev20"), FUN = list.files,
                          path = aod_dir, full.names = TRUE)
  found_files <- sapply(aer_files_dir, function(x) length(x) > 0)
  aer_files_dir = aer_files_dir[found_files]
  if(length(aer_files_dir) > 0){
    t0 <- fread(aer_files_dir[[1]], nrows = 10) # to get variable names:
    vars_aod <- intersect(grep("AOD_", names(t0), value = T), grep("nm", names(t0), value = T))
    # sort the wave lengths varnames from low to high, and update the vars_aod
    x_nm <- sort(readr::parse_number(vars_aod)) # 340, 380, ... , 1640
    vars_aod <- paste0("AOD_", x_nm,"nm") # the AOD_...nm variables
    vars_wv <- paste0("Exact_Wavelengths_of_AOD(um)_", x_nm,"nm") # the Exact wv length

    # data_path is a single file path
    read_aod <- function(data_path, date_start = NULL, date_end = NULL){
      dt <- fread(data_path, select =  c(vars0,vars_aod,vars_wv))
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
    aer_data
  } else {
    NULL
  }
}

#' Generate a tibble with a column for year and a column with a list of every
#' date in the year
#'
#' @param years vector of four-digit years
#' @return a tibble with a "year" column and a "dates" column containing a list
#'   of every date in the year
dates_year <- function(years){
  dates = lapply(years, function(y) seq.Date(as.Date(paste0(y, '-01-01')),
                                            as.Date(paste0(y, '-12-31')), 1))
  tibble::tibble(year = years, dates = dates)
}

#' Assign a month index from 1970-01-01 to all observation dates.
#'
#' Month indexes are used to batch the extraction of training data into chunks and increase targets throughput.
#' @param aer_data AERONET observation data
#' @return data.table of AERONET observation data with a monthid column added
assign_monthid <- function(aer_data){
  start = as.Date('1970-01-01')
  aer_data[, monthid := as.period(as.Date(aer_date) - start) %/% months(1)]
  aer_data
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
#' @return dataset of aer_data with pred, also with `nearest_refgrid` next to Site_Name
interpolate_aod <- function(aer_data, aer_nearest_2){
  keepcols = grep('1020|1640', names(aer_data), invert = T)
  aer_data <- aer_data[, ..keepcols]
  aer_data[,obs:=.I]
  setkey(aer_data, obs)
  vars_wv <- grep("Exact_Wavelengths_of_AOD", names(aer_data), value = T)
  vars_aod_sub <- grep("AOD_", names(aer_data), value = T)
  aer_data$N_NAAOD <- rowSums(!is.na(aer_data[,..vars_aod_sub]))
  # create long-format data
  d <- melt(aer_data[, c(vars_aod_sub, vars_wv, "obs"), with = F],
            measure = list(vars_aod_sub, vars_wv), value.name = c("aod", "exact_wv"))
  # choose non-NA and >0 AOD
  d <- d[!is.na(aod) & aod > 0]
  # ids for wavelengths observed on one side of 470nm
  obs_os <- d[, by = obs, .(oneside = (max(exact_wv)-0.47)*(min(exact_wv)-0.47)>0)][oneside ==T, obs]
  # the predict part:
  d2 <- d[, by = obs, .(AOD_470nm = exp(predict(lm(
    log(aod) ~ log(exact_wv) + I((log(exact_wv))^2)),
    newdata = data.frame(exact_wv = 0.47))))] # notice to put 0.47 not 470
  setkey(d2, obs)
  aer_data_wPred <- d2[aer_data]
  # set NA: At least 4 observations wanted for exterpolation
  aer_data_wPred[N_NAAOD<4 & obs%in%obs_os, AOD_470nm:=NA]
  # remove unncessary variables:
  aer_data_wPred <- aer_data_wPred[,-c(vars_aod_sub, vars_wv), with = F]

  # attach nearest_refgrid (also called idM21pair0 in many our script) to Aer Site_Name
  setnames(aer_data_wPred, "AERONET_Site_Name", "Site_Name")
  setkey(aer_data_wPred, Site_Name)
  aer_sites <- as.data.table(aer_nearest_2)[,c("Site_Name", "idM21pair0"), with = F]
  setkey(aer_sites, Site_Name)
  aer_data_wPred <- aer_sites[aer_data_wPred]
  setkey(aer_data_wPred, idM21pair0, stn_time)
  aer_data_wPred
}

#' Return the SF points for stations reporting data in the previously-selected
#' time period
#'
#' @param stations SF points of AERONET stations
#' @param stn_data data.table of station observations for previously-selected time period
filter_stations = function(stations, stn_data){
  stations[stations$Site_Name %in% stn_data$Site_Name, ]
}

#' Get MODIS grid unique IDs (\code{idM21pair0}) for every cell within specified
#' distance from AERONET stations.
#'
#' @param stations
#' @param refgrid_path FST with coordinates of all raster cells in the area of interest
#' @param dist_km Return cell IDs within this distance, in kilometers
cells_in_buffer = function(stations, refgrid_path, dist_km = 270){
  # set up so that function can either receive a single station (1 row) or
  # multiple stations, allowing external parallelism from dynamic branching
  gDT = read_fst(refgrid_path, as.data.table = TRUE, columns = c('idM21pair0', 'x_sinu', 'y_sinu'))

  cells_by_station = function(rownum, gDT){
    stn_coords = st_coordinates(stations[rownum, ])
    subDT = gDT[x_sinu <= stn_coords[, 'X'] + dist_km * 1000 &
                x_sinu >= stn_coords[, 'X'] - dist_km * 1000 &
                y_sinu <= stn_coords[, 'Y'] + dist_km * 1000 &
                y_sinu >= stn_coords[, 'X'] - dist_km * 1000]
    subDT[, Site_Name := stations[rownum, ]$Site_Name]
    subDT[, site_dist := as.numeric(st_distance(st_as_sf(.SD[, .(x_sinu, y_sinu)],
                                                         coords = c(1,2), crs = st_crs(stations)),
                                     stations[rownum, ])[, 1])]
    subDT[site_dist <= dist_km * 1000, .(Site_Name, idM21pair0, site_dist)]
  }
  # would exceed row maximum after ~8000 stations with 270km distance
  rbindlist(lapply(1:nrow(stations), FUN = cells_by_station, gDT = gDT))
}

#' Calculate grid X,Y positions from the top left of the most northwestern tile
#' over CONUS
#'
#' @param refgrid_path FST with coordinates of all raster cells in the area of
#'   interest
#' @return data table keyed by reference grid unique ID with X,Y offsets from
#'   northwestern corner of area of interest
calc_XY_offsets <- function(refgrid_path, ref_uid = 'idM21pair0', aoiname = 'conus'){
  if(!aoiname %in% c('conus', 'nemia')) stop('Only CONUS region has been implemented')
  refgrid = read_fst(refgrid_path, as.data.table = TRUE,
                     columns = c(ref_uid, 'tile', 'x_sinu', 'y_sinu', 'col', 'row'))

  # most northwestern tile used as origin is h08v04
  refgrid[, aoi_tile_h := as.integer(substr(tile, 2, 3)) - 8]
  refgrid[, aoi_tile_v := as.integer(substr(tile, 5, 6)) - 4]

  refgrid[, cell_x := col + 1200 * aoi_tile_h]
  refgrid[, cell_y := row + 1200 * aoi_tile_v]

  refgrid[, c('tile', 'col', 'row', 'aoi_tile_h', 'aoi_tile_v') := NULL]
}

#' Return a matrix classifying cells as being within a distance from the center.
#'
#' Assumes square matrix with odd number of rows & columns. The width (and
#' height) will be 2x the radius plus one cell.
#'
#' @param radius distance, in kilometers, from center cell to classify as
#'   inside circle
#' @param cellsize dimensions of raster cells in meters. Pixels assumed square.
#' @return a data table with offsets from the central cell that are within the
#'   specified radius
circle_mat = function(radius, cellsize = 926.6254){
  radius = radius * 1000
  width = ceiling(radius/cellsize*2+1)
  if(width%%2 == 0) width <- width + 1 # ensure odd row and column count
  height = width # assuming square matrix
  circle_df = expand.grid(mrow = 1:height, mcol = 1:width)
  center_val = median(circle_df$mrow) # assuming square matrix
  circle_df$dist_cells = sqrt((circle_df$mrow - center_val)^2 + (circle_df$mcol - center_val)^2)
  circle_df$circ = ifelse(circle_df$dist_cells * cellsize <= radius, TRUE, FALSE)
  #circle_mat = matrix(circle_df$circ, height, width)

  setDT(circle_df)
  circle_df[, offset_x := mcol - center_val]
  circle_df[, offset_y := mrow - center_val]
  circle_df[circ == TRUE, .(offset_x, offset_y)]
}

#' Roll join a single day of AERONET data to nearest MCD19A2 overpass and calculate derived
#' values from MCD19A2 within specified distances from the AERONET station.
#'
#' @param aer_data A single day of AERONET observations (provided by \code{tar_group_by})
#' @param nearby_cells
#' @param sat input "terra" or "aqua"
#' @param buffers_km vector of buffer radii in kilometers
#' @param rolldiff_limit
#' @param mcd19path path to MCD19A2 FST files
derive_mcd19_vars = function(aer_data, nearby_cells, sat,
                             buffers_km = c(10, 30, 90, 270),
                             rolldiff_limit = as.difftime(30, units = 'mins'),
                             aer_stn, mcd19path){
  setDT(aer_data)

  t1 = Sys.time()
  aer_data = interpolate_aod(aer_data, aer_stn)
  message(paste('AERONET 470nm interpolation finished in', round(Sys.time() - t1)))

  # 1. Get date & julian date
  if(aer_data[, uniqueN(aer_date)] > 1) stop('Use tar_group_by and pattern=map to send one date at a time to derive_mcd19_vars')
  this_date = aer_data[1, aer_date]
  this_daynum = format(this_date, '%j')
  this_year = year(this_date)

  # 2. Roll Join
  mcd = read_mcd19_one(sat = sat, daynum = this_daynum,
                       filepath = file.path(mcd19path, this_year),
                       load_year = this_year)
  if(!is.null(mcd)){
    # roll join will update value of the time column in X to the value of the time
    # column in i. Copy AERONET's stn_time to a column with the same name as
    # MCD19's time column so we can compare the time difference afterwards using
    # meaningful column names.
    setnames(aer_data, 'stn_time', 'overpass_time')
    aer_data[, stn_time := overpass_time]
    rj <- aer_data[mcd, roll = 'nearest', nomatch = 0] # using keys to join
    rm(aer_data)

    # Only keep AERONET to MCD19A2 joins with difference of 30 minutes or less
    rj[, rj_difftime := overpass_time - stn_time]
    rj = rj[rj_difftime <= rolldiff_limit, ]
    rj = rj[!is.na(Optical_Depth_047), ]

    # 3. Derived Values using nearby MCD19 cells
    nearby_cells = nearby_cells[Site_Name %in% rj$Site_Name]
    setnames(nearby_cells, "idM21pair0", "nearby_cellid", skip_absent=TRUE)
    # join every MODIS cell ID within 270km (default) from station
    rjbuff = rj[nearby_cells, allow.cartesian = TRUE,
                on = c(Site_Name = 'Site_Name')]

    # join mcd19 AOD values (only) to the cells in the distance buffers
    # note the join only uses idM21pair0 because this function processes one date at a time
    # if that ever changes, will need to add date to the `on` vector
    setnames(mcd, c('idM21pair0', 'Optical_Depth_047'),
                  c('nearby_cellid', 'nearby_mcd_aod'))
    mcd = mcd[, .(nearby_cellid, nearby_mcd_aod)]
    setkey(mcd, nearby_cellid)
    setkey(rjbuff, nearby_cellid)
    rjbuff_vals = mcd[rjbuff]

    calc_buff_vals <- function(distx){
      d_summary = rjbuff_vals[site_dist <= distx * 1000,
                              .(nonmissing = mean(!is.na(nearby_mcd_aod)),
                                AOD_buffer_mean = mean(nearby_mcd_aod, na.rm = T)),
                              by = idM21pair0]
      setnames(d_summary, c("nonmissing", "AOD_buffer_mean"),
               paste0(c("pNonNAAOD", "Mean_AOD"), distx, "km"))
      setkey(d_summary, idM21pair0)
    }
    new_vars = Reduce(merge, lapply(buffers_km, calc_buff_vals))
    rj_newvars = new_vars[rj]
    rm(new_vars, rj, rjbuff_vals, mcd)

    # difference between central cell MCD19 AOD and mean AOD in buffers
    for(distx in buffers_km){
      rj_newvars[, diff_AOD := Optical_Depth_047 - get(paste0('Mean_AOD', distx, 'km'))]
      setnames(rj_newvars, 'diff_AOD', paste0('diff_AOD', distx, 'km'))
    }
    if(nrow(rj_newvars) > 0){
      rj_newvars
    } else {
      # This will create an extra "V1" column in rowbound final target with all NAs,
      # but it gets around the error of writing NULL to FST format
      data.table(NA)
    }
  } else {
    data.table(NA)
  }
}

#' Open one day of MCD19A2 data from a FST file. Adds a column for the
#' satellite.
#'
#' @param sat input "terra" or "aqua", if sat !="terra", load aqua
#' @param daynum day number as character with padding to three digits, leading
#'   zeroes
#' @param filepath path to the year of MCD19A2 daily best overpasses as FST
#'   files
#' @return one day of MCD19A2 AOD data as a data.table, keyed by
#'   \code{idM21pair0, overpass_time}
read_mcd19_one <- function(sat = "terra", daynum, filepath, load_year,
                           ref_uid = 'idM21pair0',
                           columns = c("overpass_index", "Optical_Depth_047",
                                       "AOD_Uncertainty", "Column_WV", "AOD_QA",
                                       "RelAZ", "idM21pair0", "overpass_time")){
  columns = unique(c(ref_uid, 'overpass_time', columns))
  if (sat != "terra") choose = "A" else choose = "T" # load terra by default
  mcd_file = file.path(filepath, paste0('mcd19_conus_', choose, '_', load_year,
                                        daynum, '.fst'))
  if(file.exists(mcd_file)){
    dt = read.fst(mcd_file, as.data.table = TRUE, columns = columns)
    dt[, sat := choose]
    #setnames(dt, "idM21pair0", "nearest_refgrid") # rename for join later
    #setkey(dt, nearest_refgrid, join_time)
    if('overpass_time' %in% columns){
      setkeyv(dt, c(ref_uid, 'overpass_time'))
    } else {
      setkeyv(dt, ref_uid)
    }
    dt
  } else {
    NULL
  }
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
  dt[, qa_bits:= NULL]
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
  dt[, dayint:=as.integer(as.Date(overpass_time))]
  dt = dt[!is.na(get(y_var))]
  return(dt)
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
                              run_param_cv = TRUE, # run DART to select hyperparameters
                              dataXY_df = mDT[, ..features],
                              by_var = by_var,
                              n_rounds = n_rounds,
                              y_var = y_var,
                              progress = progress,
                              index_train = index_train,
                              index_test = index_test,
                              xgb_threads = xgb_threads)

  mDT <- cbind(mDT, cv_results$y_pred_dt)

  # output
  c(list(by_var = by_var, bin_list = bin_list, mDT_wPred = mDT), cv_results)
}

#' Summarize CV statistics on all initial_cv objects
#'
#' @param cv_list list of all CV output objects; may include every year and both
#'   satellites.
#' @return data.table summarizing CV statistics for each year and satellite
cv_summary <- function(cv_list){
  output = vector(mode = "list", length = length(cv_list))
  for(i in 1:length(cv_list)){
    cv_name = names(cv_list)[[i]]
    cv = cv_list[[i]]
    stats = cv_reporting(cv)
    stats$sat <- str_extract(cv_name, 'terra|aqua')
    stats$year <- substr(cv_name, 12, 15)
    output[[i]] <- stats
  }
  outDT = rbindlist(lapply(output, as.data.table))
  setkey(outDT, sat, year)
  setcolorder(outDT)
}

#' Calculate CV statistics on a single \code{initial_cv_dart()} output list object
#'
cv_reporting <- function(cv){
  dt = cv$mDT_wPred
  mae = function(v1, v2) mean(abs(v1 - v2))
  mad = function(v1) mean(abs(v1 - median(v1)))
  dt[, aod_hat := MCD19_AOD_470nm - diff_AOD_pred]

  list(
    mae_uncorrected = round(dt[, mae(MCD19_AOD_470nm, AOD_470nm)],3),
    mae_corrected   = round(dt[, mae(aod_hat, AOD_470nm)],3),
    rmse = round(cv$rmse_all_folds,3),
    mad_MCD19  = round(mad(dt$MCD19_AOD_470nm),3),
    mad_aodhat = round(mad(dt$aod_hat),3),
    stn_count = dt[, uniqueN(stn)],
    train_N = dt[, .N],
    mean_daily_overpass = round(dt[, mean(overpass_index), by = aer_date][, mean(V1)],3)
  )
}

# Prediction ####

#' Prepare prediction table
#'
#' @param pred_bbox sf bbox or named numeric vector conforming to st_bbox spec:
#'   xmin, xmax, ymax, ymin. The rectangular spatial region to predict in. Must
#'   be in MODIS sinusoidal CRS. If NULL, will operate on full CONUS.
#' @param features character vector of column names to train on
#' @param buffers_km vector of buffer radii in kilometers
#' @param refgrid_path FST with coordinates of all raster cells in the area of
#'   interest
#' @param mcd19path path to MCD19A2 FST files
#' @param aoiname the region the model was trained over ('conus' or 'nemia')
#' @param dates vector of dates to make predictions for. Must be within the same
#'   year.
#' @param sat input "terra" or "aqua"
pred_inputs <- function(pred_bbox, features, buffers_km, refgrid_path, mcd19path,
                        aoiname, sat, dates){
  if(uniqueN(year(dates)) > 1) stop('More than one year in the prediction date range')
  dates = sort(dates)

  # Get an area of interest grid including border cells needed for buffered calculations
  refgrid = calc_XY_offsets(refgrid_path, aoiname = aoiname)

  if(!is.null(pred_bbox)){
    pred_bbox = st_bbox(pred_bbox)
    m_extend = max(buffers_km) * 1000
    mod_box = st_bbox(c(xmin = -m_extend, xmax = m_extend,
                        ymin = -m_extend, ymax = m_extend))
    m_bbox = pred_bbox + mod_box
    # subset to prediction area plus an extension to support focal operations
    rgDT = refgrid[x_sinu >= m_bbox$xmin & x_sinu <= m_bbox$xmax &
                   y_sinu >= m_bbox$ymin & y_sinu <= m_bbox$ymax]
    # mark which cells to predict
    rgDT[x_sinu >= pred_bbox$xmin & x_sinu <= pred_bbox$xmax &
         y_sinu >= pred_bbox$ymin & y_sinu <= pred_bbox$ymax,
         do_preds := TRUE]
  } else {
    rgDT = refgrid
    rgDT[, do_preds := TRUE]
  }

  setkey(rgDT, cell_x, cell_y)

  # Calculate buffered values around a single cell
  buff_mcd19_vals <- function(cellid, cdf, buff_size, mcd, rgDT){
    aod_central = mcd[idM21pair0 == cellid, MCD19_AOD_470nm]
    if(length(aod_central) == 0){
      # if the requested cell is not in the FST
      return(NULL)
    } else if(is.na(aod_central)){
      # if no AOD in central cell, return a row with all NA values
      d_summary = data.table(cellid = cellid, Mean_AOD = NA, nonmissing = NA, diff_AOD = NA)
    } else {
      d_summary = tryCatch({
          buff_center = rgDT[idM21pair0 == cellid, .(cell_x, cell_y)]
          buff_offsets = cdf[, .(cell_x = offset_x + buff_center$cell_x,
                                 cell_y = offset_y + buff_center$cell_y)]
          setkey(buff_offsets)
          buff_ids = rgDT[buff_offsets, .(idM21pair0)]
          d_summary = mcd[.(buff_ids), .(cellid,
                                      'Mean_AOD' = mean(MCD19_AOD_470nm, na.rm = TRUE),
                                      'nonmissing' = sum(!is.na(MCD19_AOD_470nm))/nrow(cdf))]
          d_summary[, diff_AOD := mcd[.(cellid), MCD19_AOD_470nm] - Mean_AOD]
        },
        error = function(err){
          d_summary = data.table(cellid = cellid, Mean_AOD = NA, nonmissing = NA,
                                 diff_AOD = NA, error = as.character(err))
          setnames(d_summary, 'error', paste0('err', buff_size))
          return(d_summary)
        })
    }
    setnames(d_summary, c('Mean_AOD', 'nonmissing', 'diff_AOD', 'cellid'),
             c(paste0('Mean_AOD', buff_size, 'km'), paste0('pNonNAAOD', buff_size, 'km'),
               paste0('diff_AOD', buff_size, 'km'), 'idM21pair0'))
  }

  # Prepare MCD19A2 predictors for one day
  gather_mcd19_day <- function(this_date){
    this_year = year(this_date)
    this_daynum = format(this_date, '%j')
    # read the satellite-day MCD19A2 FST, subsetting to cells needed for prediction
    mcd = read_mcd19_one(sat = sat, daynum = this_daynum,
                         filepath = file.path(mcd19path, this_year),
                         load_year = this_year)[idM21pair0 %in% rgDT$idM21pair0]
    if(!is.null(mcd)){
      # single cell variables
      mcd[, dayint:= as.integer(as.Date(overpass_time))]
      create_qc_vars(mcd)
      setnames(mcd, "Optical_Depth_047", "MCD19_AOD_470nm")

      # buffered variables
      buff_vars = Reduce(merge,
                         # for each buffer radius:
                         lapply(buffers_km,
                                FUN = function(buff_size){
                                  # for each cell in AOI:
                                  cdf = circle_mat(buff_size)
                                  rbindlist(mclapply(rgDT[do_preds == TRUE, idM21pair0],
                                            FUN = buff_mcd19_vals,
                                            mcd = mcd, rgDT = rgDT,
                                            cdf = cdf, buff_size = buff_size,
                                            mc.cores = get.threads()), fill = TRUE)
      }))
      # join single cell variables to buffered variables
      setkey(buff_vars, idM21pair0)
      mcd[buff_vars]
    } else {
      NULL
    }
  }
  dt = rbindlist(lapply(dates, FUN = gather_mcd19_day))
  if(nrow(dt) == 0){
    dt = data.table(NA) # avoid error of writing NULL DT to FST
  } else {
    # remove unncessary columns
    dt = dt[, c('idM21pair0', features), with = FALSE]
  }
  dt
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
  model_out_path = file.path('Intermediate',
                     paste0('full_model_dart_',
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
       model_out_path = model_out_path)
}

#' Predict the difference between MCD19 and AERONET
#'
run_preds = function(data, model_file){
  data[, idM21pair0 := NULL]
  model = xgboost::xgb.load(model_file)
  predict(model, as.matrix(data))
}

#' Adjust MCD19 AOD values using predicted difference between MCD19 and AERONET.
#'
#' @return data.table of adjusted AOD values in \code{MCD19_adjust} column
adjust_mcd19 = function(data, preds){
  data[, MCD19_adjust := MCD19_AOD_470nm - preds]
}

# Maps ####

#' Compare adjusted AOD to original
#' @param refgrid_path FST with coordinates of all raster cells in the area of interest
#' @param data the data.table that was used as input to prediction function
#' @param preds numeric vector of predictions
#' @param pred_dates vector of dates that were predicted
#' @param date_index single number indicating which date to map, index in pred_dates
#' @return vertically stacked ggplots comparing original and adjusted MCD19A2 AOD
ggplot_orig_vs_adj = function(refgrid_path, data, preds, pred_dates, date_index){
  # XXX replace this date subset logic with list iteration in predinput target later
  map_dayint = as.integer(pred_dates[date_index])
  rows_per_date = nrow(data)/length(pred_dates)
  # this allows order of dates within data table to not match order of pred_dates, but that shouldn't happen
  order_dayint = which.min(abs(c(1:length(pred_dates)*rows_per_date) - max(which(data$dayint == map_dayint))))
  start_sub = rows_per_date * order_dayint - rows_per_date + 1
  data = data[start_sub:(rows_per_date * order_dayint)]
  preds = preds[start_sub:(rows_per_date * order_dayint)]

  data = adjust_mcd19(data, preds)
  rg = read_fst(refgrid_path, columns = c('idM21pair0', 'x_sinu', 'y_sinu'),
                as.data.table = TRUE)
  data[rg, c('x', 'y') := .(x_sinu, y_sinu), on = 'idM21pair0']
  orig = simple.pred.map(data, fillvar = 'MCD19_AOD_470nm')
  adj = simple.pred.map(data, fillvar = 'MCD19_adjust')
  aps = cowplot::align_plots(orig, adj, align = 'v')
  title = ggdraw() + draw_label(as.character(pred_dates[date_index]))
  cowplot::plot_grid(plotlist = list(title, aps), ncol = 1)
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
#' @param refgrid_path FST with coordinates of all raster cells in the area of interest
#' @param data the data.table that was used as input to prediction function
#' @param preds numeric vector of predictions
#' @param pred_dates
#' @param date_index
#' @param maxpixels resample raster version of predictions to approxmimately
#'   this many pixels and force the display in mapshot of this resolution
#' @return mapview object with a layer for the original and adjusted MCD19A2 AOD
mapshot_orig_vs_adj = function(refgrid_path, data, preds, pred_dates, date_index,
                               use_jenks = FALSE, maxpixels = NULL){
  # XXX replace this date subset logic with list iteration in predinput target later
  map_dayint = as.integer(pred_dates[date_index])
  rows_per_date = nrow(data)/length(pred_dates)
  # this allows order of dates within data table to not match order of pred_dates, but that shouldn't happen
  order_dayint = which.min(abs(c(1:length(pred_dates)*rows_per_date) - max(which(data$dayint == map_dayint))))
  start_sub = rows_per_date * order_dayint - rows_per_date + 1
  data = data[start_sub:(rows_per_date * order_dayint)]
  preds = preds[start_sub:(rows_per_date * order_dayint)]

  data = adjust_mcd19(data, preds)
  rg = read_fst(refgrid_path, columns = c('idM21pair0', 'x_sinu', 'y_sinu'),
                as.data.table = TRUE)
  data[rg, c('x', 'y') := .(x_sinu, y_sinu), on = 'idM21pair0']
  data = data[, .(x, y, MCD19_AOD_470nm, MCD19_adjust)]
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
  if(!dir.exists(here('Intermediate/mapshot'))) dir.create(here('Intermediate/mapshot'), recursive = T)
  mapshot_path = here('Intermediate/mapshot',
                      paste0('orig_adj_MCD19_',
                             targets:::digest_obj64(list(refgrid_path, data, preds, maxpixels)),
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
