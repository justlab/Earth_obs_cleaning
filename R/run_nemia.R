
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load Libraries                             ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(!require(here)) install.packages("here")
source(here("R/01_setup.R"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters                                 ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Output directory for rendered report
report_dir = here("reports")
if(!dir.exists(report_dir)) dir.create(report_dir)

this_year = 2017
# if you don't want to limit, do not include dates in `sel_data_bytime`
date_start = "2017-01-01" 
date_end   = "2017-12-31" 

# the aer site info. file:
aer_stn_path = "/data-coco/ECHO_PM/AeronetAODV3Level2/AOD/AOD20/aeronet_locations_v3.txt"

# the aer data files:
aer_files_path = "/data-coco/ECHO_PM/AeronetAODV3Level2/AOD/AOD20/ALL_POINTS/"

sats = c("terra", "aqua")

# MCD19 data load and save output: 
mcd19path = paste0("/data-coco/mcd19/fst/nemia/", this_year, "/")
mcd19path_CONUS = paste0("/data-coco/mcd19/fst/conus/", this_year, "/")
MCD19_T_files <- list.files(path = mcd19path_CONUS, pattern = "T", full.names = T)
MCD19_A_files <- list.files(path = mcd19path_CONUS, pattern = "A", full.names = T)
out_dir = paste0("/data-coco/mcd19/fst/conus/", this_year, "_AOD_newVar/")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# how many MCD days (files) to work on:
maxdays = as.integer(format(as.Date(paste0(this_year, "-12-31")), "%j"))
MCD_files_i = c(1:maxdays)

# for reporting purpose
y_var = "diff_AOD"
y_var_pred = "diff_AOD_pred"

# how many features to print in scatterplots
n_vars = 10

# point to cache: 
aod_cache = drake_cache(path = "/data-belle/mcd19/cache/aeronet_drake_johnathan")
set.seed(1234)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare Training Data                      ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source(here("R/drake_funcs.R"))
source(here("R/02_data_plan.R"))
(time0 <- Sys.time())
future::plan(future::multiprocess)
options(future.globals.maxSize= 100*1024^3) 
make(conus_plan, cache = aod_cache,
     parallelism = "future", jobs = 5, garbage_collection = TRUE, 
     keep_going = TRUE, memory_strategy = "preclean")
Sys.time() - time0

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run CV                                     ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source(here("R/03_cv_plan.R"))
# Allowing `add_hist` in `plot_scatter` (or use plotly) requires set lock_envir to FALSE:
system.time(
  make(conus_plan_cv, cache = aod_cache, lock_envir = FALSE)
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare Report                             ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

today = Sys.Date()
rmarkdown::render(here("AOD_CONUS_Report.Rmd"),
                  output_dir = report_dir,
                  output_file = paste0(today, "_report"),
                  output_format = rmarkdown::html_document(toc = TRUE))
