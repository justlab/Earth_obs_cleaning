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

