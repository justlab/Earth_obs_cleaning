# Functions for Aeronet's drake_plan.R 
# Part I 


# 1. Prepare Aeronet stations ------------------------------------------

#' modify_aer_stns: Modify AERONET sites (dt) to remove sites on water
#' ** Since it is a manually removal, might need to incl. more sites later **
#' @param aer_stns0 df of AERONET stations 
#' @return a new df
modify_aer_stns <- function(aer_stns0){
  aer_stns0 <- aer_stns0[!Site_Name%in%c("MVCO", "LISCO")]
  return(aer_stns0)
}

#' get_conus_buff
#' @param conus_file location of the sf CONUS file no Great Lake -- rds file
#' @return sf object of CONUS shapfile
get_conus_buff <- function(conus_file = 
  "/data-belle/LST/MODIS.LST.C6/derived/conus_GLakes_buff_sf_poly_201906.rds"){
  readRDS(conus_file)
}
# "/data-belle/LST/MODIS.LST.C6/derived/conus_sf_noGLakes.rds"

#' get_nemia_buff
#' @param states_file US states shp file
#' @return sf object of NEMIA shapfile
get_nemia_buff <- function(states_file = "/data-belle/census/states/tl_2017_us_state.shp"){
  # no longer using the old 2km buffer of NEMIA; this is same data source used for CONUS grid
  states = st_read(states_file)
  states = states[, c("STUSPS")]
  nemia = states[states$STUSPS %in% c("ME", "NH", "VT", "MA", "CT", "RI", "NY", "PA", "NJ",
                                      "DE", "MD", "DC", "VA", "WV"), ]
  nemia = st_union(nemia)
}

#' Turn Aeronet station data.table into sf project
#' @param aer_stns d.f of AERONET stations, "Site_Name", "lon", "lat", "elevm" 
#' @param sf default to TRUE, make a sf object, if FALSE, make sp object 
#' @return sf object or sp SpatialPointsDataFrame 
#' 
get_aer_spatial <- function(aer_stns, sf = TRUE){
  if(sf){
    sf::st_as_sf(aer_stns, coords = c("lon", "lat"), crs = 4326)
  }else{
    sp::SpatialPointsDataFrame(coords = aer_stns[, c("lon", "lat")], data = aer_stns, 
                           proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  }
}

#' select points by st_intersects points(aerpts) with polygon(reg_polygon)
#' @param aerpts geometry points
#' @param reg_polygon geometry polygon
#' @return aerpts: a subset of geometry points
#'
select_points <- function(aerpts, reg_polygon){
  if(st_crs(aerpts) != st_crs(reg_polygon)){ # polygons are in nad83
    aerpts <- st_transform(aerpts, crs = st_crs(reg_polygon))
  }
  aerpts <- aerpts[st_intersects(aerpts, reg_polygon, sparse = FALSE),]
  if(st_crs(aerpts) != st_crs(4326)){ aerpts <- st_transform(aerpts, crs = 4326) } # return in wgs84
  aerpts
}


#' read in reference grid, and turn into geometry points
#' @param ref_file file dir of the grid file in .fst
#' 
#' @return sf object of grid in WGS84
#' 
get_ref_grid <- function(ref_file = "/data-belle/LST/MODIS.LST.C6/derived/conus_grid_201906.fst"){
  refDT = read_fst(ref_file, as.data.table = TRUE, columns = c("idLSTpair0", "x_wgs84", "y_wgs84"))
  st_as_sf(refDT, coords = c("x_wgs84", "y_wgs84"), crs = 4326)
}



#' find points in a broad buffer, and match the nearest grid point to Aeronet sites.
#' 
#' @param conus_aer The NEMIA Aeronet sites points 
#' @param refgrid The loaded refereence grid points by \code{\function{get_ref_grid}}  
#' 
#' @return conus_aer: aer_nearest, aeronets sites with cloest grid point (idLSTpair0)
#' find near cells: find grid idLSTpair0 around 
get_nearest_cell <- function(conus_aer, refgrid){
  conus_aerNA <- st_transform(conus_aer, crs = 2163) # US National Atlas
  # limit a 1500m buffer, I agree this part is a little bit similar to `ref_in_buffer`
  aer_regions <- conus_aerNA %>% st_buffer(dist = 1500) %>% st_union
  # both the next two steps are slow:
  refgridNA <- st_transform(refgrid, crs = 2163) # US National Atlas
  int <- st_intersects(refgridNA, aer_regions, sparse = FALSE)
  refsub <- refgrid[int, ]
  
  # `refsub` is just an intermediate to speed up the process
  # to find nearby points 
  refsub <- st_transform(refsub, crs = 4326) # return in wgs84
  # match the nearest grid point to Aeronet sites
  nn_aer_ref <- st_nn(conus_aer, refsub, k = 1, returnDist = FALSE)
  # join LSTids to stations, named it as `nearest_refgrid`
  conus_aer$nearest_refgrid <- refsub[unlist(nn_aer_ref), ]$idLSTpair0
  return(conus_aer)
}

#' get nearby id names within 300 km 
#' 
#' @return geom points,df,dt, varnames: "idLSTpair0" "geometry"
get_near_cellsid <- function(conus_aer, sel_aer_region, refgrid){
  # limit the gridcell by useful aeronet sites in the first place:
  conus_aer = conus_aer[conus_aer$Site_Name%in%unique(sel_aer_region$AERONET_Site_Name),]
  conus_aerNA <- st_transform(conus_aer, crs = 2163) # US National Atlas
  aer_regions <- conus_aerNA %>% st_buffer(dist = 270000) %>% st_union
  refgridNA <- st_transform(refgrid, crs = 2163) # US National Atlas
  refsub <- refgrid[st_intersects(refgridNA, aer_regions, sparse = FALSE), ]
  return(refsub)
}


remove_site_on_water <- function(aer_nearest){
  aer_nearest[!aer_nearest$Site_Name%in%c("MVCO", "LISCO"),]
}


# 2. prepare AOD data -------------------------------------------------

# most basic variables to read
vars0 <- c("Date(dd:mm:yyyy)", "Time(hh:mm:ss)", "Day_of_Year","AERONET_Site_Name",
           "Site_Elevation(m)", "Ozone(Dobson)", "NO2(Dobson)",
           "Solar_Zenith_Angle(Degrees)", "Precipitable_Water(cm)")
# vars0_new <- c("date", "time", "dayofYear", "site",
               # "elev", "ozone", "NO2", 
               # "Solar_Zenith_Angle", "Precipitable_Water")

#' read in the Aeronet AOD measurement data from `aod20_file_dir` files.   
#' @param aod_dir where the files located on coco.   
#' @return 
#' 
get_stn_data <- function(aod_dir){
  aer_files_dir  <-  list.files(aod_dir, pattern = "*.lev20", full.names = TRUE)
  t0 <- fread(aer_files_dir[1]) # to get variable names: 
  vars_aod <- intersect(grep("AOD_", names(t0), value = T), grep("nm", names(t0), value = T))
  # sort the wave lengths varnames from low to high, and update the vars_aod
  x_nm <- sort(readr::parse_number(vars_aod)) # 340, 380, ... , 1640
  vars_aod <- paste0("AOD_", x_nm,"nm") # the AOD_...nm variables
  vars_wv <- paste0("Exact_Wavelengths_of_AOD(um)_", x_nm,"nm") # the Exact wv length
  # some other variables 
  # data_path is a single file path
  read.aod <- function(data_path){
    dt <- fread(data_path, select =  c(vars0,vars_aod,vars_wv))
    dt[, stn_time := as.POSIXct(paste(`Date(dd:mm:yyyy)`, `Time(hh:mm:ss)`), 
                           format = "%d:%m:%Y %H:%M:%S", tz = "UTC")]
    dt[, c("Date(dd:mm:yyyy)", "Time(hh:mm:ss)") := NULL]
    for (i in names(dt)) dt[get(i) == -999, (i):= NA] # set NA 
    dt
  }
  
  file_list <- lapply(aer_files_dir, read.aod) # 440s
  aer_data <- rbindlist(file_list) 
  return(aer_data)
}

#' select AOD data by station.
#' @param aer_data the aer data, limited by time, `aer_btw`
#' @param aer_stns expected to be sf points, but could be a data.table as long as it has column Site_Name
#' 
#' @return df with limited aer sites
#' 
sel_data_bystation <- function(aer_data, aer_stns){
  return(aer_data[AERONET_Site_Name %in% aer_stns$Site_Name])
}

#' select AOD data by date 
sel_data_bytime <- function(aer_data, date_start = NULL, date_end = NULL){
  if(!is.null(date_start)){ aer_data = aer_data[stn_time >= as.POSIXct(date_start, tz = "UTC"), ] }
  if(!is.null(date_end)) {aer_data = aer_data[stn_time <= as.POSIXct(date_end, tz = "UTC"),   ] }
  return(aer_data)
}


# 3. Interpolating to AOD470nm --------------------------------------------
#' prepare aer_date for interpolation of AOD470nm 
#' interpolating to wavelength 470 nm, input as 0.47 
#' @param aer_data the dataset contains all the AOD measurements and exact wavelengths
#' @return dataset of aer_data with pred, also with `nearest_refgrid` next to Site_Name
#' 
interpolate_aod <- function(aer_data, aer_nearest_2){
  aer_data <- aer_data[, -c("AOD_1020nm", "AOD_1640nm"), with = F]
  aer_data[,obs:=.I]
  setkey(aer_data, obs)
  vars_wv <- grep("Exact_Wavelengths_of_AOD", names(aer_data), value = T)
  vars_aod_sub <- grep("AOD_", names(aer_data), value = T)
  aer_data$N_NAAOD <- rowSums(!is.na(aer_data[,..vars_aod_sub]))
  # create long-format data 
  d <- melt(aer_data[, c(vars_aod_sub,vars_wv[1:22], "obs"), with = F], 
            measure = list(vars_aod_sub, vars_wv[1:22]),value.name = c("aod", "exact_wv"))
  obs_os <- d[, by = obs, .(oneside = (max(exact_wv)-0.47)*(min(exact_wv)-0.47)>0)][oneside ==T, obs]
  # choose non-NA and >0 AOD
  d <- d[!is.na(aod)]
  # the predict part: 
  d2 <- d[, by = obs, .(AOD_470nm = exp(predict(lm(log(aod) ~  log(exact_wv) + I((log(exact_wv))^2)), newdata = data.frame(exact_wv = 0.47))))] # notice to put 0.47 not 470
  setkey(d2, obs)
  aer_data_wPred <- d2[aer_data]
  # set NA: At least 4 observations wanted for exterpolation 
  aer_data_wPred[N_NAAOD<4 & obs%in%obs_os, AOD_470nm:=NA]
  # remove unncessary variables: 
  aer_data_wPred <- aer_data_wPred[,-c(vars_aod_sub, vars_wv), with = F]
  
  # attach nearest_refgrid (also called idLSTpair0 in many our script) to Aer Site_Name
  setnames(aer_data_wPred, "AERONET_Site_Name", "Site_Name")
  setkey(aer_data_wPred, Site_Name)
  aer_sites <- as.data.table(aer_nearest_2)[,c("Site_Name", "nearest_refgrid"), with = F]
  setkey(aer_sites, Site_Name)
  aer_data_wPred <- aer_sites[aer_data_wPred]
  aer_data_wPred[, join_time:=stn_time]
  setkey(aer_data_wPred, nearest_refgrid, join_time)
  aer_data_wPred
}


# 4. rolling join MCD19 (Terra, Aqua) -----------------------------------------------------------
#' read MCD19 for NEMIA
#' this function need to be replaced in the future, as for now 
#' we only read one-year data
#' @param sat input "terra" or "aqua", if sat !="terra", load aqua
read_mcd19 <- function(sat = "terra", filepath){
  if (sat != "terra") choose = "A" else choose = "T" # load terra by default
  lst_files <- list.files(path = filepath, pattern = choose, full.names = T)
  # lst_files <- list.files(path = mcd19path_CONUS, pattern = "T", full.names = T)
  ### for testing:
  lst_files <- lst_files[1:90]
  ### 
  readfile <- function(x){
    t = read.fst(x, as.data.table = T)
    # return(t[idLSTpair0%in%refsub$idLSTpair0])
    t
  }
  
  dt = rbindlist(lapply(lst_files, readfile))
  dt[,c("x", "y", "inNEMIA"):=NULL]
  dt[, join_time:=overpass_time] # join later using `overpass_time`
  dt[, sat := choose]
  setnames(dt, "idLSTpair0", "nearest_refgrid") # rename for join later 
  setkey(dt, nearest_refgrid, join_time)
  return(dt)
}

#' rolling join aer with MCD19
#' @return MCD19 rolling join Aeronet 
rolling_join <- function(mcd_sat, aer_data_wPred){
  mcd19 <- aer_data_wPred[mcd_sat, roll = 'nearest', nomatch = 0]
  # time difference
  mcd19[, diff_time_min:= as.numeric(overpass_time - stn_time)/60]
  mcd19 <- mcd19[abs(diff_time_min)<=30,]
  mcd19[, join_time:=NULL]
  return(mcd19)
}



# 5. Calculate distance ------------------------------------------------------
#' select aer sites in mcd19. further limit the stations.
#' @param region_aer all the aeronet stations
#' @param mcd19 mcd19
#' 
#' @return colocated/limited aer sites, on crs 2163
#' 
aer_in_mcd19 <- function(region_aer, mcd19){
  sel_aer = st_transform(region_aer[region_aer$Site_Name%in%unique(mcd19$Site_Name),], 
                         crs = 2163)   
  return(sel_aer)
}


#' use `st_intersects` to find grid cells in a 300km buffer.
#' This is a **slow** step, could be imporved in the future
#' 
ref_in_buffer <- function(sel_aer, refgrid){
  # not all aer sites are needed 
  radius0 <-  270000
  aod_buffers <- st_buffer(sel_aer, radius0)
  aod_buffers_list <- st_geometry(aod_buffers)
  refgrid_m <- st_transform(refgrid, crs = 2163) #sel_aer is already on crs 2163
  # # you cannot do this directly: since it will join all the buffer and produce one 
  # refsub <- refgrid_m[st_intersects(refgrid_m, aod_buffers, sparse = FALSE), ]
  join_list <- lapply(aod_buffers_list, function(x)st_intersects(refgrid_m, x))
  # list of joined I (ID) for each Aeronet station in the grid `refgrid_m`
  join_list_I <- lapply(join_list, function(x) sapply(x, function(z) if (length(z)==0) NA_integer_ else z[1]))
  # ref_list is the grid points within each circle (NA removed)
  ref_list <- lapply(join_list_I, function(x) refgrid_m[!is.na(x),])
  return(ref_list)
}

# get aer sites as a list
get_site_list <- function(sel_aer){
  site_list <- st_geometry(sel_aer)  # for calculating distance 
  # get level 1 list [i], level2 [[i]] is just a point pair
  site_list2 <- lapply(1 : length(site_list), function(i) (site_list[i])) 
  site_list2
}
# get aer site names
get_site_name <- function(sel_aer){
  n = as.list(sel_aer$Site_Name)
  n
}
# select ref.grid subset by Index for each station
select_refgrid_subset <- function(ref, p, n){
  ref$dist <- st_distance(ref, p)
  ref$Site_Name <- n
  return(as.data.table(ref))
}

#' calculate distance in buffer to each aer stn point.
#' wrap the `select_refgrid_subset` function by pmap. 
#' 
#' @return dt: a long datatable with 3 variables:
#' Site_Name, nearest_refgrid, dist_km
#' 
purrr_pmap <- function(ref, p, n){
  refDT_sub_list <- list()
  refDT_sub_list <- purrr::pmap(list(ref, p, n), select_refgrid_subset)
  # bind into one data.frame
  refDT_sub <- rbindlist(refDT_sub_list)
  refDT_sub[, dist_km:=as.numeric(dist)/1000] # change distance to km 
  # returns dt with three variables: 
  refDT_sub <- setDT(refDT_sub)[,.(Site_Name, idLSTpair0, dist_km)]
  setnames(refDT_sub, "idLSTpair0", "nearest_refgrid")
  setkey(refDT_sub, Site_Name)
  return(refDT_sub)
}


# 6. Calculate new variables ----------------------------------------------
#' Calculate new variables, For each station-day observation, we created a
#' circular buffer with a radius of 10km, 30km, 90km, and 270km around each
#' AERONET station. Then within each circular buffer, we calculated the mean MAIAC
#' AOD, the difference between MAIAC AOD and mean MAIAC AOD and the percentage of
#' non-missing MAIAC. 
#' @param aod_join_MODIS join, the rolling-joined Aeronet-MODIS
#' @param MODIS_all mcd, the original (big) MODIS dataset
#' @param refDT_sub sub, the `purrr_map` results
#' 
#' @return aod_join_MODIS with the 4x3 new variables
#' 
#' 
aod_MODIS_newVars <- function(aod_join_MODIS, MODIS_all, refDT_sub){
  aod_join_MODIS <- aod_join_MODIS[!is.na(Optical_Depth_047),]
  aod_join_MODIS[, day:=format(stn_time, "%Y-%m-%d")]
  setkey(aod_join_MODIS, Site_Name)
  # join by site_name 
  setnames(refDT_sub, "nearest_refgrid", "buf_refgrid", skip_absent=TRUE)
  aod_join_buffer <- aod_join_MODIS[refDT_sub, allow.cartesian = T]
  message("The cartesian joined dataset between collocated Aeronet with MODIS sites data and base grid in each circles, dim is: ", 
        paste(dim(aod_join_buffer), collapse = " x "))
  
  # join back to the MODIS file using `buf_refgrid` in the joined dataset
  MODIS_all[, day:=format(join_time, "%Y-%m-%d")]
  setnames(MODIS_all, c("nearest_refgrid", "Optical_Depth_047"), 
           c("buf_refgrid", "buffer_AOD_470"), skip_absent = T)
  setkey(aod_join_buffer, buf_refgrid, day)
  setkey(MODIS_all, buf_refgrid, day)
  # buffer_AOD_470 is the AOD in the buffer grid cells
  # aod_join_buffer2 get AOD for all the buffer grid cells
  aod_join_buffer2 <- MODIS_all[,.(buf_refgrid, day, buffer_AOD_470)][aod_join_buffer]
  
  # calculate variables for 4 different distances 
  # for each nearest_refgrid (which associates with Site_Name) in the joined dataset
  dist0 <- c(10, 30, 90, 270)
  calculate.vars.in.buffer <- function(distx){
    aod_join_buffer_new <- aod_join_buffer2[dist_km < distx, 
                                            .(nonmissing = mean(!is.na(buffer_AOD_470)), 
                                              AOD_buffer_mean = mean(buffer_AOD_470, na.rm = T)),
                                            by =.(nearest_refgrid, day)]
    setnames(aod_join_buffer_new, c("nonmissing", "AOD_buffer_mean"), 
             paste0(c("pNonNAAOD", "Mean_AOD"),distx, "km"))
    setkey(aod_join_buffer_new, nearest_refgrid, day)
    return(aod_join_buffer_new)
  }
  new_varlist <- invisible(lapply(dist0, calculate.vars.in.buffer))
  # joined by column 
  new_vars <- Reduce(merge, new_varlist)
  setkey(new_vars, nearest_refgrid, day)
  setkey(aod_join_MODIS, nearest_refgrid, day)
  aod_join_MODIS <- new_vars[aod_join_MODIS]
  
  # calculate `diff_AOD`: AOD - mean AOD
  calculate.diff <- function(distx){
    aod_join_MODIS[,diff_AOD:= get(paste0("Mean_AOD",distx,"km")) - AOD_Uncertainty]
    setnames(aod_join_MODIS, "diff_AOD", paste0("diff_AOD",distx,"km"))
  }
  invisible(lapply(dist0, calculate.diff))
  return(aod_join_MODIS)
}


# 7. CV -------------------------------------------------------------------
# to create qc variables 
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

#' run the cross-validation by station with RFE.
#' 
#' This function is the 4th layer wrap of Kodi's random search cv function
#' I apologize for the confusion if anyone is planning to dig into it. 
#' Please feel free to contact me. 
#' 
#' 
run_cv <- function(dt){
  setnames(dt, "Optical_Depth_047", "MCD19_AOD_470nm", skip_absent=TRUE)
  # The dependent variable: diff_AOD = MCD19 - AERONET = Optical_Depth_047 - AOD_470nm
  dt[, diff_AOD := MCD19_AOD_470nm - AOD_470nm]
  
  # The predictors
  features = c("MCD19_AOD_470nm",
               "dayint", 
               "AOD_Uncertainty", "Column_WV", "RelAZ",
               "qa_best", 
               do.call(paste0,expand.grid(
                 c("pNonNAAOD", "Mean_AOD", "diff_AOD"),
                                          paste0(c(10, 30, 90, 270),"km"))))
  

  dt <- create_qc_vars(dt)
  dt[, dayint:=as.integer(as.Date(day))]
  
  # cv
  cv_results <- run.k.fold.cv.rfe.wrap(
                modeldt1 = dt, 
                stn_var = "Site_Name",
                features0 = features, 
                sat = "", y_var = "diff_AOD", 
                run_param_cv = T, run_rfe = T)
  return(cv_results)
}


# 8. Report ---------------------------------------------------------------
report_rmse <- function(rferesults){
  report.r.squared(dataXY = rferesults$modeldt1_wPred)
}

plot_rmse <- function(rferesults){
  rfe.rmse.plot(rferesults$rmse_rfe)
}

plot_AOD <- function(rferesults){
  dt = rferesults$modeldt1_wPred
  p1 = scatter.plot.diagonal(dt, x = "MCD19_AOD_470nm", y = "AOD_470nm")
  p2 = scatter.plot.diagonal(dt, x = "MCD19_AOD_470nm", y = "diff_AOD")
  p3 = scatter.plot.diagonal(dt, x = "AOD_470nm", y = "diff_AOD")
  p4 = scatter.plot.diagonal(dt, x = "MCD19_AOD_470nm", y = "diff_AOD_pred")
  list(p1,p2,p3,p4)
}

plot_scatter <- function(rferesults){
  dt = rferesults$modeldt1_wPred
  var_list <- rferesults$features_rank_rfe[1:n_vars] #(features)
  fig_list <- lapply(var_list, scatter.plot.simple, y = "diff_AOD", data = dt)
  fig_list
}

get_SHAP_long <- function(rferesults){
  var_list <- rferesults$features_rank_rfe #(features)
  # plot SHAP
  shap_long <- shap.prep(shap_contrib = rferesults$shap_score, 
                         X_train = rferesults$modeldt1_wPred[,..var_list])
  return(shap_long)
}

plot_SHAP <- function(shap_long){
  shap.plot.summary(shap_long, scientific = T)
}

plot_SHAP_scatter <- function(rferesults, shap_long){
  var_list <- rferesults$features_rank_rfe[1:n_vars]
  fig_list <- lapply(var_list, function(x)shap.plot.dependence.color(shap_long, x = x, y = x, 
                                                                     color_feature = "MCD19_AOD_470nm"))
  fig_list
}
