# Functions for aeronet_drake.R ~~~~~~~~~~~~~~~~~~~~~~ ####


# 1. Prepare subset of Aeronet stations ------------------------------------------

#' modify_aer_stns: Modify AERONET sites (dt) to remove sites on water
#' @param aer_stns0 df of AERONET stations 
#' @return a new df
modify_aer_stns <- function(aer_stns0){
  aer_stns0 <- aer_stns0[!Site_Name%in%c("MVCO", "LISCO")]
  return(aer_stns0)
}

#' get_conus_buff
#' @param conus_file location of the sf CONUS file no Great Lake -- rds file
#' @return sf object of CONUS shapfile
get_conus_buff <- function(conus_file = "/data-belle/LST/MODIS.LST.C6/derived/conus_sf_noGLakes.rds"){
  readRDS(conus_file)
}

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
#' @return sf object of grid in WGS84
#' 
get_ref_grid <- function(ref_file = "/data-belle/LST/MODIS.LST.C6/derived/conus_grid_201906.fst"){
  refDT = read_fst(ref_file, as.data.table = TRUE, columns = c("idLSTpair0", "x_wgs84", "y_wgs84"))
  st_as_sf(refDT, coords = c("x_wgs84", "y_wgs84"), crs = 4326)
}

#' read in the Aeronet AOD measurement data from `data_path`
get_stn_data <- function(data_path, keep_columns){
  dt = fread(data_path, select = keep_columns, fill = TRUE, skip = 6)
  dt[, day := as.POSIXct(paste(`Date(dd:mm:yyyy)`, `Time(hh:mm:ss)`), 
                         format = "%d:%m:%Y %H:%M:%S", tz = "UTC")]
  dt[, c("Date(dd:mm:yyyy)", "Time(hh:mm:ss)") := NULL]
  return(dt)
}


#' find points in buffer within dist
#' @param aerpts The NEMIA Aeronet sites points 
#' @param refpts The loaded referecen grid points by \code{\function{get_ref_grid}}  
#' @return candidate_refpts
#' 
points_in_buffer <- function(aerpts, refpts, dist = 1500){
  aerptsNA <- st_transform(aerpts, crs = 2163) # US National Atlas
  aer_regions <- aerptsNA %>% st_buffer(dist = dist) %>% st_union
  refptsNA <- st_transform(refpts, crs = 2163) # US National Atlas
  refsub <- refpts[st_intersects(refptsNA, aer_regions, sparse = FALSE), ]
  st_transform(refsub, crs = 4326) # return in wgs84
}


#' match the nearest grid point to Aeronet sites  
#' @param aerpts The NEMIA Aeronet sites points 
#' @param candidate_refpts all the points_in_buffer
get_nearest_cell <- function(aerpts, candidate_refpts){
  nn_aer_ref <- st_nn(aerpts, candidate_refpts, k = 1, returnDist = FALSE)
  # join LSTids to stations
  aerpts$nearest_refgrid <- candidate_refpts[unlist(nn_aer_ref), ]$idLSTpair0
  aerpts
}
remove_site_on_water <- function(aer_nearest){
  aer_nearest[!aer_nearest$Site_Name%in%c("MVCO", "LISCO"),]
}


# 2. prepare AOD measurement data -------------------------------------------------
# select AOD data by station and date 
sel_data_bystation <- function(aer_data, aer_stns){
  # aer_stns expected to be sf points, but could be a data.table as long as it has column Site_Name
  aer_data[AERONET_Site %in% aer_stns$Site_Name]
}


sel_data_bytime <- function(aer_data, date_start = NULL, date_end = NULL){
  if(!is.null(date_start)){ aer_data = aer_data[day >= as.POSIXct(date_start, tz = "UTC"), ] }
  if(!is.null(date_end)) {aer_data = aer_data[day <= as.POSIXct(date_end, tz = "UTC"),   ] }
  return(aer_data)
}
