# Functions for aeronet_drake.R ~~~~~~~~~~~~~~~~~~~~~~ ####

get_stn_data <- function(data_path, keep_columns){
  dt = fread(data_path, select = keep_columns, fill = TRUE, skip = 6)
  dt[, day := as.POSIXct(paste(`Date(dd:mm:yyyy)`, `Time(hh:mm:ss)`), 
                         format = "%d:%m:%Y %H:%M:%S", tz = "UTC")]
  dt[, c("Date(dd:mm:yyyy)", "Time(hh:mm:ss)") := NULL]
  return(dt)
}
sel_data_bytime <- function(aer_data, date_start = NULL, date_end = NULL){
  if(!is.null(date_start)){ aer_data = aer_data[day >= as.POSIXct(date_start, tz = "UTC"), ] }
  if(!is.null(date_end)){   aer_data = aer_data[day <= as.POSIXct(date_end, tz = "UTC"),   ] }
  return(aer_data)
}
get_ref_grid <- function(ref_file = "/data-belle/LST/MODIS.LST.C6/derived/conus_grid_201906.fst"){
  refDT = read_fst(ref_file, as.data.table = TRUE, columns = c("idLSTpair0", "x_wgs84", "y_wgs84"))
  st_as_sf(refDT, coords = c("x_wgs84", "y_wgs84"), crs = 4326)
}
get_conus_buff <- function(conus_file = "/data-belle/LST/MODIS.LST.C6/derived/conus_sf_noGLakes.rds"){
  readRDS(conus_file)
}
get_nemia_buff <- function(states_file = "/data-belle/census/states/tl_2017_us_state.shp"){
  # no longer using the old 2km buffer of NEMIA; this is same data source used for CONUS grid
  states = st_read(states_file)
  states = states[, c("STUSPS")]
  nemia = states[states$STUSPS %in% c("ME", "NH", "VT", "MA", "CT", "RI", "NY", "PA", "NJ",
                                      "DE", "MD", "DC", "VA", "WV"), ]
  nemia = st_union(nemia)
}
get_aer_spatial <- function(aer_stns, sf = TRUE){
  if(sf){
    st_as_sf(aer_stns, coords = c("lon", "lat"), crs = 4326)
  }else{
    SpatialPointsDataFrame(coords = aer_stns[, c("lon", "lat")], data = aer_stns, 
                           proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  }
}
select_points <- function(aerpts, reg_polygon){
  if(st_crs(aerpts) != st_crs(reg_polygon)){ # polygons are in nad83
    aerpts <- st_transform(aerpts, crs = st_crs(reg_polygon))
  }
  aerpts <- aerpts[st_intersects(aerpts, reg_polygon, sparse = FALSE),]
  if(st_crs(aerpts) != st_crs(4326)){ aerpts <- st_transform(aerpts, crs = 4326) } # return in wgs84
  aerpts
}
points_in_buffer <- function(aerpts, refpts, dist = 1500){
  aerptsNA <- st_transform(aerpts, crs = 2163) # US National Atlas
  aer_regions <- aerptsNA %>% st_buffer(dist = dist) %>% st_union
  refptsNA <- st_transform(refpts, crs = 2163) # US National Atlas
  refsub <- refpts[st_intersects(refptsNA, aer_regions, sparse = FALSE), ]
  #st_transform(refsub, crs = 4326) # return in wgs84
}
get_nearest_cell <- function(aerpts, candidate_refpts){
  nn_aer_ref <- st_nn(aerpts, candidate_refpts, k = 1, returnDist = FALSE)
  # join LSTids to stations
  aerpts$nearest_refgrid <- candidate_refpts[unlist(nn_aer_ref), ]$idLSTpair0
  aerpts
}
sel_data_bystation <- function(aer_data, aer_stns){
  # aer_stns expected to be sf points, but coule be a data.table as long as it has column Site_Name
  aer_data[AERONET_Site %in% aer_stns$Site_Name]
}

# mcd19 processing ####

# open matching years of mcd19 data, Aqua and Terra if 2002+, extracting data for the closest cells to AERONET stations


# join_mcd19(mcd19path, date_start, date_end, aer_nearest, aer_data){
#   # get years in the date range
#   
#   # loop over years (xxx append year to path)
#   
#   # do Terra
#     # matching files, get matching cells for all
#     # do the join - check times
#   
#   # if 2002+, do Aqua (adding columns to Terra table)
#     # do in parallel with Terra? Make this loop an expanded set of targets by year?
#   
#    
# }
# get_mcd19_day(path, cells){
#   
# }
