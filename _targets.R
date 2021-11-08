library(targets)
library(tarchetypes)
library(future)
# library(future.callr)
# plan(callr)
#plan(multicore)
plan(multisession)
options(xgb.threads = 8) # used by xgboost_cv_RFE

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Configuration ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tar_option_set(
  packages = c('lubridate',
               'stringr',
               'log4r',
               'data.table',
               'fst',
               'sf',
               'magrittr',
               'future',
               'here',
               'raster',
               'ggplot2',
               'SHAPforxgboost',
               'Just.universal',
               'xgboost',
               'tibble'),
  format = 'qs',
  workspace_on_error = TRUE)

tar_config_set(store = '/data-belle/cache/aod_targets/')

source('R/globals.R')
source('R/drake_funcs.R')
source('R/functions.R')
source('R/xgboost_cv_RFE.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Targets ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

process_years = 2010:2019
region_values = list(regions = aoiname)
date_table = dates_year(process_years)
sat_values = list(sat = sats)
buffers_km = c(10, 30, 90, 270)
features = c("MCD19_AOD_470nm", "dayint", "AOD_Uncertainty",
             "Column_WV", "RelAZ", "qa_best",
             do.call(paste0, expand.grid(
               c("pNonNAAOD", "Mean_AOD", "diff_AOD"),
               paste0(buffers_km, "km"))))

set1_targets = list(
  tar_target(mcd19path,
             '/data-coco/mcd19/fst/conus_full'),
  tar_target(aer_stn_path,
             '/data-coco/ECHO_PM/AeronetAODV3Level2/AOD/AOD20/aeronet_locations_v3.txt',
             format = 'file'),
  tar_target(aer_files_path,
             '/data-coco/ECHO_PM/AeronetAODV3Level2/AOD/AOD20/ALL_POINTS/',
             format = 'file'),
  tar_target(refgrid_path,
             '/data-belle/LST/MYD21A1/derived/conus_grid_2020.fst',
             format = 'file'),
  tar_target(refras_path,
             '/data-belle/LST/MYD21A1/derived/conus_myd21_stack.tif',
             format = 'file'),
  tar_target(aer_stations,
             fread(aer_stn_path, col.names = c('Site_Name', 'lon', 'lat', 'elevm')),
             format = 'fst_dt'),

  # Area of interest ####
  tar_map( # region mapping
    values = region_values,
    tar_target(buff,
               get_aoi_buffer(regions)),
    tar_target(aer,
               select_stations(aer_stations, buff, refgrid_path, refras_path)),
    tar_target(nearby_cells,
               cells_in_buffer(aer, refgrid_path),
               format = 'fst_dt'),

    # Load AERONET data ####
    tar_target(all_dates,
               lubridate::as_date(unlist(date_table$dates))),
    tar_target(aer_nospace,
               sf::st_drop_geometry(aer)),
    tar_group_by(aer_bystation,
                 aer_nospace, Site_Name),
    tar_target(aer_data,
               get_stn_data(aod_dir = aer_files_path, stations = aer_bystation),
               pattern = map(aer_bystation)),
    tar_target(aer_filtered,
               filter_aer_bydate(aer_data, all_dates),
               format = 'fst_dt'),
    tar_target(aer_months,
               assign_monthid(aer_filtered),
               format = 'fst_dt'),
    tar_group_by(aer_bymonth,
                 aer_months, monthid),

    # Load MCD19A2 AOD ####
    tar_map( # sat mapping
      values = sat_values,
      tar_target(mcd19_vars,
                 purrr::map_dfr(aer_bymonth %>% split(.$aer_date),
                                derive_mcd19_vars,
                                buffers_km = buffers_km,
                                nearby_cells = nearby_cells,
                                sat = sat,
                                aer_stn = aer_nospace,
                                mcd19path = mcd19path),
                 pattern = map(aer_bymonth),
                 format = 'fst_dt',
                 storage = 'worker'),

      # Model ####
      tar_map( # date range (year) mapping
        unlist = FALSE,
        values = date_table,
        names = 'year',
        tar_target(modelinput,
                   prepare_dt(mcd19_vars, date_range = dates),
                   format = 'fst_dt'),
        tar_target(initial_cv,
                   initial_cv_dart(modelinput,
                     y_var = "diff_AOD",
                     features = features,
                     stn_var = "Site_Name"))
      ),
      # later, move within year mapping
      # Prediction ####
      tar_target(pred_dates, as.Date('2018-06-01')),
      tar_target(predinput, pred_inputs(
        pred_bbox = dc_sinu,
        features = features,
        buffers_km = buffers_km,
        refgrid_path = refgrid_path,
        mcd19path = mcd19path,
        aoiname = regions,
        sat = sat,
        dates = pred_dates),
        pattern = map(pred_dates),
        format = 'fst_dt',
        storage = 'worker')
    )
  )
)

# Combine all initial CV results into a list ####
set1u = unlist(set1_targets)
combined_target = tar_combine(combined_cv, set1u[grep('initial_cv', names(set1u))],
                              command = list(!!!.x))

# Render CV report ####
cv_report_target = tar_render(initial_cv_report, 'R/initial_cv_report.Rmd',
                       params = list(cv_names = names(combined_cv)))


# Final targets list ####
list(set1_targets, combined_target, cv_report_target)
