
#' Return AERONET points that intersect the polygon describing the area of
#' interest, and include the cell ID (idM21pair0) from the reference CONUS MODIS
#' grid.
#' @param stations data.table including coordiantes (lon, lat) points of AERONET stations
#' @param reg_polygon SF polygon defining the area of interest
#' @return a subset of geometry points within the given polygon
select_stations <- function(stations, reg_polygon){
  aerpts = st_as_sf(stations, coords = c("lon", "lat"), crs = 4326)
  if(st_crs(aerpts) != st_crs(reg_polygon)){ # polygons are in nad83
    aerpts <- st_transform(aerpts, crs = st_crs(reg_polygon))
  }
  aerpts <- aerpts[st_intersects(aerpts, reg_polygon, sparse = FALSE),]
  aerpts <- station_cell_ids(aerpts) # now in sinusoidal CRS
  aerpts
}

#' Get the unique ID of the grid cell an AERONET station is in.
#' GLOBALS: refgrid_path, refras_path
#' @param stations_sf SF points of AERONET stations
#' @return SF points of AERONET stations with unique MODIS grid cell IDs added (idM21pair0)
station_cell_ids = function(stations_sf){
  gras = raster(refras_path, band = 4)
  stations_sf = st_transform(stations_sf, crs = st_crs(gras))
  stations_sf$cell_index = raster::extract(gras, as(stations_sf, 'Spatial'))
  
  gDT = read_fst(refgrid_path, as.data.table = TRUE, columns = c('idM21pair0', 'cell_index'))
  setkey(gDT, cell_index)
  setDT(stations_sf)
  setkey(stations_sf, cell_index)
  stations_sf[gDT, idM21pair0 := idM21pair0]
  stations_sf = st_sf(stations_sf)
}

