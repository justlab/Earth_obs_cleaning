# See `drake_funcs.R` for functions used in drake plan #

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

# Libraries
suppressPackageStartupMessages({
  library(sf)
  library(magrittr)
  library(drake)
  library(data.table)
  library(fst)
  library(future)
  library(here)
  #install.packages("/home/rushj03/dev/nngeo", repos = NULL, type = "source")
  # library(nngeo, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.6")
  #library(sp)
  #library(rgdal)
  #library(rnaturalearth) # REPLACE WITH CONUS FILE USED IN GRID XXX
  #library(rgeos)
  library(nngeo)
  #library(nngeo, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.6")
})

# Parameters ####

aer_stn_path = "/data-coco/ECHO_PM/AeronetAODV3Level2/AOD/AOD20/aeronet_locations_v3.txt"

aer_data_path = '~/qnap_geo/aeronet/All_Sites_Times_All_Points_AOD20_20180603.dat'
aer_data_path = '/scratch/temp/All_Sites_Times_All_Points_AOD20_20180603.dat' # don't have time to play with Ethernet all day
# the aer data with AOD470nm
aer_data_path_new = "~/Intermediate/AOD_data/190628_aod20_wPred_NEMIA_idLST.fst"
  
aer_files_path = "/data-coco/ECHO_PM/AeronetAODV3Level2/AOD/AOD20/ALL_POINTS/"
date_start = "2018-01-01" # do not include in sel_data_bytime if don't want to limit
date_end   = "2018-12-31" # do not include in sel_data_bytime if don't want to limit
geo_region = "NEMIA"
sats = c("terra", "aqua")
# MCD19 data: 
mcd19path = "/data-coco/mcd19/fst/nemia" # now there is only 2018 data 
# point to cache: 
ssd_cache = new_cache(path = "/scratch/cache/aeronet_drake") 

# Drake Plans  ------------------------------------------------------------
# load functions for drake: 
source(here("R/drake_funcs.R"))
# new plan Yang Liu
nemia_plan <- drake_plan(
  # 1.  Aeroent ---------------------------------------------------------------------
  
  # AERONET station locations file: 
  aer_stns = fread(file = file_in(aer_stn_path), col.names = c("Site_Name", "lon", "lat", "elevm")),
  conus_buff = get_conus_buff(),
  nemia_buff = get_nemia_buff(),
  aerpts = get_aer_spatial(aer_stns),               # wgs84
  conus_aer = select_points(aerpts, conus_buff),    # AERONET sites in CONUS as points
  nemia_aer = select_points(conus_aer, nemia_buff), # AERONET sites in NEMIA as points
  
  # get nearest MODIS cells 
  refgrid = get_ref_grid(),
  candidate_refpts = points_in_buffer(nemia_aer, refgrid),
  aer_nearest = get_nearest_cell(nemia_aer, candidate_refpts), # Aeronet matched to grid, 174x4
  aer_nearest_NEMIA = remove_site_on_water(aer_nearest), # remove 2 sites 
  
  # measurements. Editing values of date_start and date_end should make use of the cache
  # if range was previously calc'd
  aer_data = get_stn_data(aod_dir = file_in(aer_files_path)),  # 26171206x56, 11G
  aer_btw = target(sel_data_bytime(aer_data, date_start, date_end),
                   transform = map(date_start = !!date_start, date_end = !!date_end)),
  
  # selected measurements, transform by same date range if not automatic in name
  nemia = target(sel_data_bystation(aer_btw, nemia_aer),
                      transform = map(aer_btw)), # use map to map to aer_btw
  
  # 2. Interpolation -----------------------------------------------------------
  # so instead of running all the Aeronet data, only run what is used
  # it also join with the nearest grid id 
  wPred = target(interpolate_aod(nemia, aer_nearest_NEMIA),transform = map(nemia))
)
nemia_plan
#drake_plan_source(nemia_plan)
nemia_config <- drake_config(nemia_plan, cache = ssd_cache) # show the dependency
# vis_drake_graph(nemia_config, from = names(nemia_config$layout))
vis_drake_graph(nemia_config)

# outdated(nemia_config)
future::plan("multicore")
make(nemia_plan, cache = ssd_cache) # write update into the cache 
cached(cache = ssd_cache)
loadd(cache = ssd_cache) # load all object 

# to see the updated version
nemia_config <- drake_config(nemia_plan, cache = ssd_cache) # show the dependency
vis_drake_graph(nemia_config)


# Join MCD19 --------------------------------------------------------------
aer_data_wPred <- readd("wPred_nemia_aer_btw_2018.01.01_2018.12.31", cache = ssd_cache) # load all object 
plan1 <- drake_plan(
  # 3. WPred Join MCD19 -----------------------------------------------------------
  mcd = target(read_lst(sat), transform = map(sat = c("terra", "aqua"))),
  join = target(rolling_join(mcd_sat = mcd, aer_data_wPred),
                transform = map(mcd))
)
config2 <- drake_config(plan1, cache = ssd_cache) # show the dependency
vis_drake_graph(config2)
time0 <- Sys.time()
make(plan1, cache = ssd_cache) # write update into the cache 
Sys.time() - time0
cached(cache = ssd_cache)
loadd(cache = ssd_cache) # load all object 
join_mcd_aqua = readd("join_mcd_aqua", cache = ssd_cache)


