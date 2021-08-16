library(targets)
library(tarchetypes)
library(future)
# plan(multisession)
library(future.callr)
plan(callr)

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
               'Just.universal',
               'tibble'),
  format = 'qs')

tar_resources_fst(compress = 100)
tar_config_set(store = '/data-belle/cache/aod_targets/')

source('R/globals.R')
source('R/drake_funcs.R')
source('R/functions.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Targets ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

region_values = list(regions = aoiname)
time_values = list(date_start = date_start, date_end = date_end)
sat_values = list(sat = sats)

list(
  tar_target(aer_stations, 
             fread(aer_stn_path, col.names = c("Site_Name", "lon", "lat", "elevm")),
             format = 'fst_dt'),
  
  tar_map( # region mapping
    values = region_values, 
    
    tar_target(buff, get_aoi_buffer(regions)),
    tar_target(aer, select_stations(aer_stations, buff)),
    tar_target(nearby_cells, cells_in_buffer(aer),
               format = 'fst_dt'),
    
    tar_map( # time mapping
      values = time_values,
      tar_target(aer_data, get_stn_data(aod_dir = aer_files_path, stations = aer,
                                        date_start, date_end),
                 format = 'fst_dt'),
      tar_target(aer_filtered, filter_stations(aer, aer_data)),
      tar_group_by(aer_bydate, aer_data, aer_date),
      tar_map( # sat mapping
        values = sat_values,
        tar_target(mcd19_vars, derive_mcd19_vars(aer_bydate, nearby_cells, sat), 
                   pattern = head(map(aer_bydate), n = 20),
                   format = 'fst_dt',
                   storage = 'worker')
      )
    )
  )
)
