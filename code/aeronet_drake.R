# beginning of script where we will have a plan to:

# select Aeronet in CONUS
# select Aeronet in NEMIA
# calculate nearest idLSTpair0's. May need shortcut that uses buffer around stations to reduce possibility set. 

# select from the daily data the subsets above (one of the two)
# calc columns as in old script
# as.numeric columns as necessary
# remove rows without AOD? probably not - they may have CWV
# ? which other columns to we want?
# restrict to 2000+ (parameter so that different subsets could be taken later)

# Parameters ~~~~~~~~~~~~~~~~~~~~~~ ####

aer_stn_path = '~/qnap_geo/aeronet/aeronet_locations_v3_20180604.txt'
aer_data_path = '~/qnap_geo/aeronet/All_Sites_Times_All_Points_AOD20_20180603.dat'
date_start = "2000-01-01" # set to NULL if don't want to limit
date_end   = "2018-12-31" # set to NULL if don't want to limit
geo_region = "NEMIA"
aer_columns = c("AERONET_Site", "Date(dd:mm:yyyy)", "Time(hh:mm:ss)", "Day_of_Year", "AOD_440nm", 
                "Site_Latitude(Degrees)", "Site_Longitude(Degrees)", # remove later, curious if there are mobile stations
                "Solar_Zenith_Angle(Degrees)", "Precipitable_Water(cm)")

ssd_cache = new_cache(path = "/scratch/cache/aeronet_drake") # currently on Belle

# Functions  ~~~~~~~~~~~~~~~~~~~~~~ ####

get_stn_data <- function(data_path, date_start, date_end, keep_columns){
  dt = fread(data_path, select = keep_columns)
  dt[, day := as.POSIXct(paste(`Date(dd:mm:yyyy)`, `Time(hh:mm:ss)`), 
                                 format = "%d:%m:%Y %H:%M:%S", tz = "UTC")]
  dt[, c("Date(dd:mm:yyyy)", "Time(hh:mm:ss)") := NULL]
  if(!is.null(date_start)){ dt = dt[day >= as.Date(date_start), ] }
  if(!is.null(date_end)){   dt = dt[day <= as.Date(date_end),   ] }
  return(dt)
}

# Drake Plans  ~~~~~~~~~~~~~~~~~~~~ ####

nemia_plan <- drake_plan(
  aer_stns = fread(file = file_in(aer_stn_path), col.names = c("Site_Name", "lon", "lat", "elevm")),
  aer_data = target(get_stn_data(file_in(aer_data_path), date_start, date_end, aer_columns),
             transform = map(date_start = !!date_start, date_end = !!date_end)),
  
)

nemia_plan
drake_plan_source(nemia_plan)

nemia_config <- drake_config(nemia_plan, cache = ssd_cache)
vis_drake_graph(nemia_config, from = names(nemia_config$layout))

outdated(nemia_config)
make(nemia_plan, cache = ssd_cache)
cached(cache = ssd_cache)




