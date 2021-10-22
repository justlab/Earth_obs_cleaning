
#' Return AERONET points that intersect the polygon describing the area of
#' interest, and include the cell ID (idM21pair0) from the reference CONUS MODIS
#' grid.
#' @param stations data.table including coordiantes (lon, lat) points of AERONET stations
#' @param reg_polygon SF polygon defining the area of interest
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
#' @param stations_sf SF points of AERONET stations
#' @return SF points of AERONET stations with unique MODIS grid cell IDs added (idM21pair0)
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
#' `tar.gz` file. Only load the data for the specified stations.
#' @param aod_dir directory where the `.lev20` files were extracted to.
#' @param stations SF points of AERONET stations
#' @param date_start subset to observations this date or later
#' @param date_end subset to observations this date or earlier
#' @return data.table of AERONET observations, joined with MODIS reference grid unique ID
#'
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


#' Generate a tibble with a column for year and a column with a list of every date in the year
#'
#' @param years vector of four-digit years
#' @return a tibble with a "year" column and a "dates" column containing a list of every date in the year
#'
dates_year <- function(years){
  dates = lapply(years, function(y) seq.Date(as.Date(paste0(y, '-01-01')),
                                            as.Date(paste0(y, '-12-31')), 1))
  tibble::tibble(year = years, dates = dates)
}

#' Assign a month index from 1970-01-01 to all observation dates.
#'
#' Month indexes are used to batch the extraction of training data into chunks and increase targets throughput.
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


#' prepare aer_date for interpolation of AOD470nm
#' interpolating to wavelength 470 nm, input as 0.47
#' @param aer_data the dataset contains all the AOD measurements and exact wavelengths
#' @return dataset of aer_data with pred, also with `nearest_refgrid` next to Site_Name
#'
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
#' @param stations SF points of AERONET stations
#' @param stn_data data.table of station observations for previously-selected time period
filter_stations = function(stations, stn_data){
  stations[stations$Site_Name %in% stn_data$Site_Name, ]
}

#' Get MODIS grid unique IDs (idM21pair0) for every cell within specified
#' distance from AERONET stations.
#' @param stations
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

#' Roll join a single day of AERONET data to nearest MCD19A2 overpass and calculate derived
#' values from MCD19A2 within specified distances from the AERONET station.
#' @param aer_data A single day of AERONET observations (provided by tar_group_by)
#' @param nearby_cells
#' @param sat
#' @param buffers_km
#' @param rolldiff_limit
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
    aer_data[, overpass_time := stn_time]
    rj <- aer_data[mcd, roll = 'nearest', nomatch = 0,
                   on = c(idM21pair0 = 'idM21pair0', overpass_time = 'overpass_time')]
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
    setnames(mcd, 'idM21pair0', 'nearby_cellid')
    rjbuff_vals = mcd[, .(nearby_cellid = nearby_cellid,
                          nearby_mcd_aod = Optical_Depth_047)][
                       rjbuff, on = c(nearby_cellid = 'nearby_cellid')]

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

#' Open one day of MCD19A2 data from a FST file. Adds a column for the satellite.
#' @param sat input "terra" or "aqua", if sat !="terra", load aqua
#' @param daynum day number as character with padding to three digits, leading zeroes
#' @param filepath path to the year of MCD19A2 daily best overpasses as FST files
read_mcd19_one <- function(sat = "terra", daynum, filepath, load_year){

  if (sat != "terra") choose = "A" else choose = "T" # load terra by default
  mcd_file = file.path(filepath, paste0('mcd19_conus_', choose, '_', load_year,
                                        daynum, '.fst'))
  if(file.exists(mcd_file)){
    dt = read.fst(mcd_file, as.data.table = TRUE,
                  columns = c("overpass_index", "Optical_Depth_047",
                              "AOD_Uncertainty", "Column_WV", "AOD_QA", "RelAZ",
                              "idM21pair0", "overpass_time"))

    dt[, sat := choose]
    #setnames(dt, "idM21pair0", "nearest_refgrid") # rename for join later
    #setkey(dt, nearest_refgrid, join_time)
    setkey(dt, idM21pair0, overpass_time)
    dt
  } else {
    NULL
  }
}


# TEMPORARY BARELY-CHANGED CV PREP FUNCTIONS FROM DRAKE
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


# Do initial CV with hyperparameter selection with DART and rank features
# Depends on functions in xgboost_cv_RFE.R
initial_cv_dart <- function(
  data,
  #sat = '',
  y_var,
  features,
  k_fold = 5,
  n_rounds = 100,
  stn_var = NULL,
  day_var = NULL,
  #run_param_cv = TRUE,
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
  # removed: will use hyperparams as a parameter to RFE function
  # # write param list so in rfe process there is no need to run param search
  # xgb_param_list_full <- cv_results$xgb_param_list2
  # if (!dir.exists(here("Intermediate"))) dir.create(here("Intermediate"))
  # saveRDS(xgb_param_list_full, here("Intermediate", "xgb_param_list_full.rds"))

  mDT <- cbind(mDT, cv_results$y_pred_dt)

  # output
  c(list(by_var = by_var, bin_list = bin_list, mDT_wPred = mDT), cv_results)
}

cv_reporting <- function(cv){
  dt = cv$mDT_wPred
  mae = function(v1, v2) mean(abs(v1 - v2))
  mad = function(v1) mean(abs(v1 - median(v1)))
  dt[, aod_hat := MCD19_AOD_470nm - diff_AOD_pred]

  list(
    mae_uncorrected = dt[, mae(MCD19_AOD_470nm, AOD_470nm)],
    mae_corrected   = dt[, mae(aod_hat, AOD_470nm)],
    mad_MCD19  = mad(dt$MCD19_AOD_470nm),
    mad_aodhat = mad(dt$aod_hat)
  )
}
