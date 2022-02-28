# Functions not currently in use
# Johnathan finds it useful to have these available for consulting without
# having to search through git history, but feel free to remove this file after
# he's departed.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# removed from pred_inputs()
# These were used to calculate focal statistics using data.table
# It had been tuned and profiled, but was still slower than terra focal stats, at least for large areas of interest

refgrid = calc_XY_offsets(refgrid_path, aoiname = aoiname)
# then handle the pred_bbox for labeling do_preds in rgDT as currently done...
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# removed from gather_mcd19_day(), replaced with call to buff_mcd19_raster()

# Old data.table focal stats:
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lookups between original and aggregated resolutions now done with Terra

#' Return a data.table that assoicates cells of the original AOD resolution and
#' the aggregated resolution.
#'
#' Contains columns with the reference cell ID from the pair function (provided
#' as `ref_uid` parameter), cell ID from the original reference raster
#' (`cell_index`), cell ID from the newly created raster (`cell_index_new`), the
#' row and column of the raster (`cell_row`, `cell_col`), and columns for the
#' aggregated version of the raster (`agg_cell_index`, `agg_cell_row`,
#' `agg_cell_col`).
#'
make_agg_lookup <- function(mcd_refras, agg_level, refgrid_path, ref_uid = 'idM21pair0'){
  # create and stack reference rasters
  sr_refras = rast(mcd_refras)
  agg_cellids <- terra::aggregate(sr_refras, fact = agg_level)
  agg_cellids[] <- 1:ncell(agg_cellids)
  names(agg_cellids) <- 'agg_cell_index'
  agg_cellids[['agg_cell_row']] <- rowFromCell(agg_cellids, 1:ncell(agg_cellids))
  agg_cellids[['agg_cell_col']] <- colFromCell(agg_cellids, 1:ncell(agg_cellids))
  disagg_cellids = terra::disagg(agg_cellids, fact = agg_level) %>%
    terra::crop(., sr_refras)
  rm(agg_cellids)
  sr_refras[['cell_index_new']] <- 1:ncell(sr_refras)
  sr_refras[['cell_row']] <- rowFromCell(sr_refras, 1:ncell(sr_refras))
  sr_refras[['cell_col']] <- colFromCell(sr_refras, 1:ncell(sr_refras))
  stacked = rast(list(sr_refras, disagg_cellids))
  rm(sr_refras, disagg_cellids)

  # convert to data.table and join ref_uid
  ref_agg = setDT(as.data.frame(stacked))
  rg = read_fst(refgrid_path, as.data.table = TRUE,
                columns = c(ref_uid, 'cell_index'))
  ref_agg[rg, c(ref_uid) := get(ref_uid), on = 'cell_index']
  ref_agg
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# No longer opening mcd19 FSTs output from MAIACprocessing project. Now opening
# data from HDFs via VRTs.

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# No longer looking up MCD19 data from data tables

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
