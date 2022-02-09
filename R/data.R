# load HDFs of MCD19A2

# FOR TESTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(targets)
library(data.table)
library(terra)
library(stringr)
library(jsonlite)
fvec = list.files('~/qnap_geo/MCD19A2/HDF/2010.01.31', pattern = '.hdf$', full.names = TRUE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# extract orbit times from HDFs
orbit_info <- function(hdfpath){
  md = fromJSON(describe(hdfpath, options = 'json'))
  ocount = as.numeric(md$metadata[[1]]$Orbit_amount)
  otimes = str_split(md$metadata[[1]]$Orbit_time_stamp, '  ')[[1]][1:ocount]
  list(tile = rep(str_extract(hdfpath, 'h[:digit:][:digit:]v[:digit:][:digit:]'), ocount),
       layer_index = 1:ocount,
       sat = str_sub(otimes, -1),
       orbit_time = as.POSIXct(str_sub(otimes, 1, nchar(otimes)-1), tz = 'UTC', format = '%Y%j%H%M'),
       hdf = hdfpath)
}

# make a DT of orbits binned by time for every tile in a day
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

