# load HDFs of MCD19A2

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
  oDT[sat == 'A', sat := 'aqua']
  oDT[sat == 'T', sat := 'terra']
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

#' Get numeric indices for requested HDF subdatasets
get_sd_ids <- function(hdf_path){
  md = fromJSON(describe(hdf_path, options = 'json'))
  md = md$metadata$SUBDATASETS
  md = md[grep('DESC', names(md))]
  ids = as.integer(str_extract(names(md), '(?<=SUBDATASET_).*(?=_DESC)'))
  sd_names = str_extract(md, '(?<=] ).*(?= g)')
  data.table(id = ids, sdname = sd_names)
}

#' Create a list of SpatRasters from a list of paths to VRTs.
#'
#' @param op_vrts is the list of VRTs by subdatset from `get_overpasses_vrts()`,
#'   subset to this overpass
#' @param opDT data.table from `bin_overpasses()`, subset to this overpass
#' @param with_time create a mosaic of single-cell rasters representing the
#'   overpass time per HDF tile.
#' @return a list with the same length as there are unique raster resolutions
#'   among the subdatasets. For MCD19A2, it will be a list of length three when
#'   `with_time = TRUE`. The first item contains 4 rasterized VRTs for the high
#'   resolution subdatasets, then 1 raster with lower res, and finally 1 time
#'   raster with a single cell per HDF tile. These terra SpatRaster objects can
#'   not be serialized as-is, since they are pointers. They must be consumed by
#'   the same process.
rasterize_vrts <- function(op_vrts, opDT, with_time = TRUE){

  # rasterize the subdatasets and group by resolution
  ras_by_res = list()
  vrt_rasters = lapply(op_vrts, rast)
  resvec = sapply(vrt_rasters, function(x) res(x)[1])
  for(this_res in unique(resvec)){
    ras_by_res[[paste0('r', floor(this_res))]] <-
      rast(vrt_rasters[which(resvec == this_res)])
  }

  if(with_time){
    # Create time rasters by tile
    time_rasters = list()
    for(t in 1:nrow(opDT)){
      r1 = rast(opDT[t, hdf], subds = 1, lyrs = 1)
      time_rasters[[t]] <- rast(ext(r1), crs = crs(r1), nrows = 1, ncols = 1,
                                vals = as.numeric(opDT[t, orbit_time]),
                                names = 'overpass_time')
    }
    time_merged = terra::merge(sprc(time_rasters))
    ras_by_res[['time']] <- time_merged
  }
  ras_by_res
}

#' Get all overpass mosaics as VRTs by date.
#'
#' @param hdf_paths character vector of paths to all HDF tiels for the date.
#' @param binDT data.table from `bin_overpasses()`
#' @param load_sat 'aqua' or 'terra'. Default 'aqua'.
#' @param vrt_path directory to store overpass VRTs that reference HDF files.
#' @param load_sds character vector of which subdatasets to load. All layers
#'   (overpasses) will be loaded for each of the subdatasets. The default
#'   subdatasets are set to the those being used for MCD19A2.
#' @return paths to VRTs in a nested list by overpass and subdataset for this date
get_overpasses_vrts <- function(hdf_paths, binDT, load_sat = sats, vrt_path,
  load_sds = c('Optical_Depth_047', 'AOD_Uncertainty', 'Column_WV', 'AOD_QA', 'RelAZ')){
  load_sat = match.arg(load_sat)
  binDT = binDT[sat == load_sat]

  sdids = get_sd_ids(hdf_paths[1]) # assuming identical for all HDFs
  sdids = sdids[sdname %in% load_sds]

  # list of all overpasses
  op_vrts = list()
  for(op_id in unique(binDT$overpass_bin)){
    # list of subdatasets in this overpass
    sd_vrts = list()
    for(sd_row in 1:nrow(sdids)){
      sd_vrts[[sd_row]] <- custom_vrt_day(load_sat, sdids[sd_row, id],
                                          sdids[sd_row, sdname], binDT,
                                          op_id, vrt_path)
    }
    names(sd_vrts) <- sdids$sdname
    op_vrts[[paste0('op_bin_', op_id)]] <- sd_vrts
  }
  op_vrts
}

#' Create a directory to store overpass VRTs
prepare_vrt_directory <- function(base_path){
  vrt_path = file.path(base_path, 'vrts')
  if(!dir.exists(vrt_path)) dir.create(vrt_path)
  vrt_path
}

# lyrs being a vector from the layer_index column from bin_overpasses

#' Create VRT files for a single overpass and single subdataset that references
#' the appropriate layer in multiple HDF files.
#'
#' @details tested with system GDAL 3.0.4 released 2020-01-28
#' @param this_sat 'aqua' or 'terra'.
#' @param sdnum numeric index of the HDF subdataset.
#' @param sdname name of the HDF subdataset.
#' @param binDT data.table from `bin_overpasses()`.
#' @param op_bin overpass bin number.
#' @param vrt_path directory to store overpass VRTs that reference HDF files.
custom_vrt_day <- function(this_sat, sdnum, sdname, binDT, op_bin, vrt_path){
  binDT = binDT[overpass_bin == op_bin & sat == this_sat, ]
  oldwd = getwd()
  hdf_path = dirname(binDT$hdf[1])
  setwd(vrt_path)
  this_date = as.Date(binDT[1, orbit_time])
  vrtname = paste(this_sat, this_date, sdname, op_bin, 'vrt', sep = '.')
  # write VRT with all SourceBands to 1
  system2('gdalbuildvrt', paste('-b 1', '-sd', sdnum, '-overwrite', vrtname,
                                paste(binDT$hdf, collapse = ' ')),
          stdout = NULL)

  # edit each tile in the output VRT's SourceBands to match specified overpass
  newbands = paste0('      <SourceBand>', binDT$layer_index, '</SourceBand>')
  vrt = readLines(vrtname)
  matches = grep('      <SourceBand>1</SourceBand>', vrt)
  for(i in 1:length(matches)){
    vrt[matches[i]] <- newbands[i]
  }
  writeLines(vrt, vrtname)
  setwd(oldwd)
  file.path(vrt_path, vrtname)
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
