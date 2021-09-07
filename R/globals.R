# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Globals ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Region
aoiname = "conus"

# Reference grid
refgrid_path = "/data-belle/LST/MYD21A1/derived/conus_grid_2020.fst"
refras_path = "/data-belle/LST/MYD21A1/derived/conus_myd21_stack.tif"

# Year
this_year = 2017
# if you don't want to limit, do not include dates in `sel_data_bytime`
date_start = paste0(this_year, "-01-01")
#date_end = paste0(this_year, "-01-02") # for testing
date_end = paste0(this_year, "-12-31") 

# the aer site info. file:
aer_stn_path = "/data-coco/ECHO_PM/AeronetAODV3Level2/AOD/AOD20/aeronet_locations_v3.txt"

# the aer data files:
aer_files_path = "/data-coco/ECHO_PM/AeronetAODV3Level2/AOD/AOD20/ALL_POINTS/"

# which MODIS platforms to use
sats = c("terra", "aqua")

# MCD19 data load and save output: 
mcd19path_CONUS = paste0("/data-coco/mcd19/fst/conus_full/", this_year, "/")
# out_dir = paste0("/data-coco/mcd19/fst/conus/", this_year, "_AOD_newVar/")
# if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# for reporting purpose
y_var = "diff_AOD"
y_var_pred = "diff_AOD_pred"

# how many features to print in scatterplots
n_vars = 10

set.seed(1234)
