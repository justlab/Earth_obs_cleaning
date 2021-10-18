# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Globals ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Region
aoiname = "conus"

# which MODIS platforms to use
sats = c("terra", "aqua")

# out_dir = paste0("/data-coco/mcd19/fst/conus/", this_year, "_AOD_newVar/")
# if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# for reporting purpose
y_var = "diff_AOD"
y_var_pred = "diff_AOD_pred"

# how many features to print in scatterplots
n_vars = 10

set.seed(1234)
