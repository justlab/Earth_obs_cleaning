# load HDFs of MCD19A2

# FOR TESTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(targets)
library(data.table)
library(terra)
library(stringr)
library(jsonlite)
library(mapview)
library(sf)
#fvec = list.files('~/qnap_geo/MCD19A2/HDF/2010.01.31', pattern = '.hdf$', full.names = TRUE)
hdf_root = '/home/rushj03/qnap_geo/MCD19A2/HDF'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Orbit Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#' Extract orbit times from HDFs
orbit_info <- function(hdf_path){
  md = fromJSON(describe(hdf_path, options = 'json'))
  ocount = as.numeric(md$metadata[[1]]$Orbit_amount)
  otimes = str_split(md$metadata[[1]]$Orbit_time_stamp, '  ')[[1]][1:ocount]
  list(tile = rep(str_extract(hdf_path, 'h[:digit:][:digit:]v[:digit:][:digit:]'), ocount),
       layer_index = 1:ocount,
       sat = str_sub(otimes, -1),
       orbit_time = as.POSIXct(str_sub(otimes, 1, nchar(otimes)-1), tz = 'UTC', format = '%Y%j%H%M'),
       hdf = hdf_path)
}

#' Make a DT of orbits binned by time for every tile in a day
bin_overpasses <- function(hdf_paths){
  oDT = rbindlist(lapply(hdf_paths, FUN = orbit_info))
  setkey(oDT, sat, orbit_time)
  oDT[, time_lag := data.table::shift(.SD, 1), by = sat, .SDcols = 'orbit_time']
  oDT[, time_diff := orbit_time - time_lag]
  oDT[, time_lag := NULL]
  for(this_sat in unique(oDT$sat)){
    cuts = oDT[sat == this_sat & (is.na(time_diff) | time_diff > as.difftime(1800, units = 'secs')), orbit_time]
    oDT[sat == this_sat, overpass_bin := findInterval(orbit_time, cuts)]
  }
  oDT[, .(tile, sat, layer_index, overpass_bin, orbit_time, hdf)]
}

# Data Loading ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#' Stack the requested subdataset layers, disaggregating coarser resolution
#' rasters to the finest resolution.
#'
#' @param hdf_path full path to a single HDF of satellite data.
#' @param load_sds character vector of which subdatasets to load.
#' @param load_layers numeric vector of indices of layers to load within
#'   subdatasets (which correspond to desired overpasses).
#' @param time_layer_value POSIXct time. Optional. If provided, will add a numeric layer
#'   to the raster containing the time as seconds since 1970.
stack_and_disagg <- function(hdf_path, load_sds, load_layers, time_layer_value = NULL){
  bysds = lapply(load_sds, FUN = rast, x = hdf_path, lyrs = load_layers)
  # Disaggregate the coarser layer (RelAZ in MCD19A2) to the finer layers' resolution.
  # This method assumes square pixels.
  resvec = sapply(bysds, function(x) res(x)[1])
  disagg_level = round(resvec / min(resvec))
  for(i in which(disagg_level > 1)){
    bysds[[i]] <- disagg(bysds[[i]], disagg_level[i])
  }
  names(bysds) <- load_sds
  if(!is.null(time_layer_value)){
    if(!'POSIXct' %in% class(time_layer_value)) stop('Provide the time as a POSIXct object')
    newindex = length(bysds)+1
    bysds[[newindex]] <- init(bysds[[1]], as.numeric(time_layer_value))
    names(bysds[[newindex]]) <- 'overpass_time'
  }
  rast(bysds)
}

#' Get all overpass mosaics by date.
#'
#' File organization is expected to match USGS server organization (see
#' parameters). This function does not return wrapped/packed `terra` rasters, so
#' should not be used for a target. It should be used within another target
#' where the output remains in memory.
#' @param hdf_root Path to directory containing date subdirectories formatted
#'   YYYY.MM.DD. The date subdirectories contain the satellite HDFs.
#' @param load_date Date object for date to load.
#' @param load_sat 'A' for Aqua or 'T' for Terra. Default Aqua.
#' @param load_sds character vector of which subdatasets to load. All layers
#'   (overpasses) will be loaded for each of the subdatasets. The default
#'   subdatasets are set to the those being used for MCD19A2.
get_daily_overpasses <- function(hdf_root, load_date, load_sat = c('A', 'T'),
  load_sds = c('Optical_Depth_047', 'AOD_Uncertainty', 'Column_WV', 'AOD_QA', 'RelAZ')){
  load_sat = match.arg(load_sat)

  hdf_paths = list.files(file.path(hdf_root, format(load_date, '%Y.%m.%d')),
                         pattern = '\\.hdf$', full.names = TRUE)
  dt = bin_overpasses(hdf_paths)
  dt = dt[sat == load_sat]

  # outer loop by overpass_bin
  merge_overpass = function(bin){
    obdt = dt[overpass_bin == bin]
    # loop by tile
    tlist = vector(mode = 'list', length = nrow(obdt))
    for(i in 1:nrow(obdt)){
      trow = obdt[i]
      tlist[[i]] <- stack_and_disagg(trow$hdf, load_sds, trow$layer_index,
                                     time_layer_value = trow$orbit_time)
    }
    obin <- terra::merge(sprc(tlist))
  }
  olist = lapply(unique(dt$overpass_bin), merge_overpass)
  # olist = mclapply(unique(dt$overpass_bin), merge_overpass, mc.cores = 6)
  olist
}

#' Convert a single SpatRaster (stack of all layers) to data.table
overpass_to_table <- function(sr){
  dt = setDT(as.data.frame(sr, xy = TRUE, na.rm = FALSE))
  setnames(dt, c('x', 'y'), c('x_sinu', 'y_sinu'))
  dt
}

# Mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#' Make a raster with overpass bins for a day (used for map only)
hdf_layer_as_binary <- function(hdf_path, layer_index, overpass_bin,
                                sd_name = c('AOD_QA', 'Optical_Depth_047')){
  sd_name = match.arg(sd_name)
  r1 = rast(hdf_path, subds = sd_name, lyrs = layer_index)
  set.values(r1, overpass_bin, cells = which(!is.na(r1[])))
  r1
}

#' Get densified polygonal HDF tile borders
tile_boundaries <- function(hdf_paths){
  st_envelope = function(x) st_as_sfc(st_bbox(x))
  ras_names = str_extract(hdf_paths, 'h[:digit:][:digit:]v[:digit:][:digit:]')
  rasters <- lapply(hdf_paths, rast, subds = 'AOD_QA')
  envelopes <- sapply(rasters, st_envelope)
  env_sf <- st_sf(st_sfc(envelopes))
  st_crs(env_sf) = st_crs(rasters[[1]])
  # Create more points along the rectangular bounding box so that the borders
  # can be reprojected with a curve. Densification param is in map units.
  dens_env = densify(vect(env_sf), 100000)
  sf_dens = st_sf(tiles = ras_names, st_geometry(st_as_sf(dens_env)))
  sf_dens
}

#' Create an interactive mapview of MCD19A2 overpasses in a day
#' @param hdf_paths character vector of full paths for all tiles in a day to map
#' @param map_sat 'A' for Aqua or 'T' for Terra. Default Aqua.
mapview_overpasses <- function(hdf_paths, map_sat = c('A', 'T')){
  map_sat = match.arg(map_sat)
  sat_name = if(map_sat == 'T') 'Terra' else 'Aqua'

  # Get mosaic of overpasses
  dt = bin_overpasses(hdf_paths)
  dt = dt[sat == map_sat]

  olist = vector(mode = 'list', length = length(unique(dt$overpass_bin)))
  for(bin in unique(dt$overpass_bin)){
    obdt = dt[overpass_bin == bin]
    tlist = vector(mode = 'list', length = nrow(obdt))
    for(i in 1:nrow(obdt)){
      trow = obdt[i]
      tlist[[i]] <- hdf_layer_as_binary(trow$hdf, trow$layer_index, trow$overpass_bin,
                                         sd_name = 'AOD_QA')
    }
    obin <- terra::merge(sprc(tlist))
    olist[[bin]] <- obin
  }
  # Overpasses are sequential east to west, where they overlap, the value will
  # be the mean of the two overpasses. This method works for CONUS, but breaks
  # down closer to the poles where there are more overlapping overpasses!
  omos = terra::mosaic(sprc(olist), fun = 'mean')

  # Tile boundaries to overlay on raster
  tb = tile_boundaries(hdf_paths)

  # qualitative 12 class from ColorBrewer2, in dark,light order instead of light,dark
  raster_colors = c("#1f78b4", "#a6cee3", "#33a02c", "#b2df8a", "#e31a1c", "#fb9a99",
                    "#ff7f00", "#fdbf6f", "#6a3d9a", "#cab2d6", '#b15928', '#ffff99')
  # put in temporary data.table to allow recycling for >12 overpasses
  raster_vals = unique(omos)
  suppressWarnings({
    sel_colors = data.table(raster_vals, raster_colors)[1:nrow(raster_vals), raster_colors]
  })
  mapview(raster::raster(omos), layer.name = paste(sat_name, 'Overpasses'),
          na.color = '#BEBEBE00', method = 'ngb', col.regions = sel_colors) +
    mapview(tb, alpha.regions = 0, color = 'black', legend = FALSE)
}
