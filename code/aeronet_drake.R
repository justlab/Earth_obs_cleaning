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
aer_data_path = '/scratch/temp/All_Sites_Times_All_Points_AOD20_20180603.dat' # don't have time to play with Ethernet all day
date_start = "2000-01-01" # do not include in sel_data_bytime if don't want to limit
date_end   = "2018-12-31" # do not include in sel_data_bytime if don't want to limit
geo_region = "NEMIA"
aer_columns = c("AERONET_Site", "Date(dd:mm:yyyy)", "Time(hh:mm:ss)", "Day_of_Year", "AOD_440nm", 
                "Site_Latitude(Degrees)", "Site_Longitude(Degrees)", # remove later, curious if there are mobile stations
                "Solar_Zenith_Angle(Degrees)", "Precipitable_Water(cm)")

ssd_cache = new_cache(path = "/scratch/cache/aeronet_drake") # currently on Belle

# Functions  ~~~~~~~~~~~~~~~~~~~~~~ ####

get_stn_data <- function(data_path, keep_columns){
  dt = fread(data_path, select = keep_columns, fill = TRUE, skip = 6)
  dt[, day := as.POSIXct(paste(`Date(dd:mm:yyyy)`, `Time(hh:mm:ss)`), 
                                 format = "%d:%m:%Y %H:%M:%S", tz = "UTC")]
  dt[, c("Date(dd:mm:yyyy)", "Time(hh:mm:ss)") := NULL]
  return(dt)
}
sel_data_bytime <- function(aer_data, date_start = NULL, date_end = NULL){
  if(!is.null(date_start)){ aer_data = aer_data[day >= as.POSIXct(date_start, tz = "UTC"), ] }
  if(!is.null(date_end)){   aer_data = aer_data[day <= as.POSIXct(date_end, tz = "UTC"),   ] }
  return(aer_data)
}

# Drake Plans  ~~~~~~~~~~~~~~~~~~~~ ####

nemia_plan <- drake_plan(
  aer_stns = fread(file = file_in(aer_stn_path), col.names = c("Site_Name", "lon", "lat", "elevm")),
  aer_data = get_stn_data(file_in(aer_data_path), aer_columns),
  aer_btw = target(sel_data_bytime(aer_data, date_start, date_end),
                    transform = map(date_start = !!date_start, date_end = !!date_end))
)

# note: line 25564581 out of 26171185 ends early. possibly others, but that's where fread warned
# 329 s for aer_data last time from QNAP, 108s from SSD

nemia_plan
drake_plan_source(nemia_plan)

nemia_config <- drake_config(nemia_plan, cache = ssd_cache)
vis_drake_graph(nemia_config, from = names(nemia_config$layout))
vis_drake_graph(nemia_config)

outdated(nemia_config)
make(nemia_plan, cache = ssd_cache)
cached(cache = ssd_cache)




