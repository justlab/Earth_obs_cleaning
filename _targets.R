library(targets)
library(tarchetypes)
tar_option_set(
  packages = c('data.table',
               'fst',
               'sf',
               'magrittr',
               'future',
               'here',
               'raster',
               'ggplot2',
               'SHAPforxgboost',
               'Just.universal'),
  format = 'qs')
#tar_option_set(debug = 'buff')
tar_resources_fst(compress = 100)
tar_config_set(store = '/data-belle/cache/aod_targets/')

source('R/globals.R')
source('R/drake_funcs.R')
source('R/functions.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Targets ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

values = list(regions = aoiname)

tar_map(
  values = values, 
  
  tar_target(aer_stations, 
             fread(aer_stn_path, col.names = c("Site_Name", "lon", "lat", "elevm")),
             format = 'fst_dt'),
  tar_target(buff, get_aoi_buffer(regions)),
  tar_target(aer, select_stations(aer_stations, buff))
)
