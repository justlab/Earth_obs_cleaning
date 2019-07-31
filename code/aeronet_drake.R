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

# Libraries  ~~~~~~~~~~~~~~~~~~~~~~ ####
library(sf)
library(magrittr)
library(drake)
library(data.table)
library(fst)
library(future)
#install.packages("/home/rushj03/dev/nngeo", repos = NULL, type = "source")
library(nngeo)
#library(nngeo, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.6")

# Parameters ~~~~~~~~~~~~~~~~~~~~~~ ####

aer_stn_path = '~/qnap_geo/aeronet/aeronet_locations_v3_20180604.txt'
aer_data_path = '~/qnap_geo/aeronet/All_Sites_Times_All_Points_AOD20_20180603.dat'
aer_data_path = '/scratch/temp/All_Sites_Times_All_Points_AOD20_20180603.dat' # don't have time to play with Ethernet all day
#date_start = "2000-01-01" # do not include in sel_data_bytime if don't want to limit
date_start = "2018-01-01" # do not include in sel_data_bytime if don't want to limit
date_end   = "2018-12-31" # do not include in sel_data_bytime if don't want to limit
geo_region = "NEMIA"
aer_columns = c("AERONET_Site", "Date(dd:mm:yyyy)", "Time(hh:mm:ss)", "Day_of_Year", "AOD_440nm", 
                #"Site_Latitude(Degrees)", "Site_Longitude(Degrees)", # confirmed stations do not move!
                "Solar_Zenith_Angle(Degrees)", "Precipitable_Water(cm)")

mcd19path = "/data-coco/mcd19/fst/nemia" 
ssd_cache = new_cache(path = "/scratch/cache/aeronet_drake") # currently on Belle

# Functions ~~~~~~~~~~~~~~~~~~~~ 

source("code/aeronet_processing.R")

# Drake Plans  ~~~~~~~~~~~~~~~~~~~~ ####

nemia_plan <- drake_plan(
  # station locations
  aer_stns = fread(file = file_in(aer_stn_path), col.names = c("Site_Name", "lon", "lat", "elevm")),
  conus_buff = get_conus_buff(),
  nemia_buff = get_nemia_buff(),
  aerpts = get_aer_spatial(aer_stns),               # wgs84
  conus_aer = select_points(aerpts, conus_buff),    
  nemia_aer = select_points(conus_aer, nemia_buff), 
  
  # get nearest MODIS cells 
  refgrid = get_ref_grid(),
  candidate_refpts = points_in_buffer(nemia_aer, refgrid),
  aer_nearest = get_nearest_cell(nemia_aer, candidate_refpts),
  
  # measurements. Editing values of date_start and date_end should make use of the cache if range was previously calc'd
  aer_data = get_stn_data(file_in(aer_data_path), aer_columns),
  aer_btw = target(sel_data_bytime(aer_data, date_start, date_end),
                   transform = map(date_start = !!date_start, date_end = !!date_end)),
  
  # selected measurements, transform by same date range if not automatic in name
  nemia_data = target(sel_data_bystation(aer_btw, nemia_aer),
                      transform = map(aer_btw)),
  
  # joined mcd19 MAIAC values. transform by same date range if not automatic in name
  
)



nemia_plan
#drake_plan_source(nemia_plan)
nemia_config <- drake_config(nemia_plan, cache = ssd_cache)
#vis_drake_graph(nemia_config, from = names(nemia_config$layout))
vis_drake_graph(nemia_config)

outdated(nemia_config)
plan("multicore")
make(nemia_plan, cache = ssd_cache)
cached(cache = ssd_cache)

# testing ####
loadd(candidate_refpts, cache = ssd_cache)
loadd(nemia_aer, cache = ssd_cache)
mapview(candidate_refpts, col.regions = "blue") + mapview(nemia_aer, col.regions = "red")
loadd(aer_nearest, cache = ssd_cache)
nemia_data = readd("nemia_data_aer_btw_2018.01.01_2018.12.31", cache = ssd_cache)

nemia_data
nemia_data[, uniqueN(AERONET_Site)] # 19
nemia_data[, tdiff := day - shift(day), by = AERONET_Site] 
nemia_tdiff = nemia_data[, .(obscount = .N, mean_tdiff = mean(tdiff, na.rm = T), 
              med_tdiff = median(tdiff, na.rm = T)), by = AERONET_Site]
library(ggplot2)
ggplot(nemia_tdiff) + geom_density(aes(mean_tdiff))
ggplot(nemia_tdiff[mean_tdiff<1000]) + geom_density(aes(mean_tdiff))


# notes ####
# 2000-2018: 174 unique NEMIA sites, 130 in the data

# old notes ####

## doesn't look like separate plans are possible; maybe fork and edit drake_plan() sometime
# data_plan <- drake_plan(
#   nemia_plan, # would need to expand this to the rows
#   aer_data = get_stn_data(file_in(aer_data_path), aer_columns),
#   aer_btw = target(sel_data_bytime(aer_data, date_start, date_end),
#                    transform = map(date_start = !!date_start, date_end = !!date_end)) )


# note: line 25564581 out of 26171185 ends early. possibly others, but that's where fread warned
# 329 s for aer_data last time from QNAP, 108s from SSD
