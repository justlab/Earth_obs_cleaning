Sys.setenv(RENV_CONFIG_STARTUP_QUIET = "TRUE")
renv::load()

suppressPackageStartupMessages(
   {library(data.table)
    library(targets)})

config = yaml::read_yaml("config.yaml")

data.dir = config$data.dir
stopifnot(dir.exists(data.dir))
geonexl2.dir = file.path(data.dir, "geonexl2")
intermediate.path = function(...)
   file.path(data.dir, 'intermediate', ...)
dir.create(intermediate.path(), showWarnings = F)
download = function(from, to, ...)
    download.update.meta(from, file.path(data.dir, "downloads"), to, ...)
satellite_hdf_root = file.path(data.dir, 'earthdata')
dir.create(satellite_hdf_root, showWarnings = F)

n.workers = config$n.workers
if (is.null(n.workers))
    n.workers = max(1, parallel::detectCores() - 2)

# The most general workflow-defining variables.
Wf = as.environment(config$workflow)
stopifnot(Wf$outcome %in% c("aod"))
stopifnot(
    Wf$satellite.product == "mcd19a2" &&
       Wf$satellite %in% c("terra", "aqua") ||
    Wf$satellite.product == "geonexl2" &&
       Wf$satellite == "goes16")
stopifnot(Wf$ground.product %in% c("aeronet"))
stopifnot(
    # The region can be a special string…
    Wf$region == "conus" ||
    # …or a filepath…
    startsWith(Wf$region, "/") ||
    # …or a URL, ending with an appropriate file extension.
    startsWith(Wf$region, "http:") ||
    startsWith(Wf$region, "https:"))

# True configurability of other `Wf` items is not implemented now
# but might be added later.

Wf$years = switch(Wf$satellite.product,
    mcd19a2 = 2000 : 2021,
    geonexl2 = 2018 : 2019)
Wf$dates = seq(
    lubridate::make_date(min(Wf$years)),
    lubridate::make_date(max(Wf$years), 12, 31),
    by = 1)
Wf$date.example = switch(Wf$satellite.product,
  # This should be a date for which the satellite data of interest
  # exists on all tiles.
    mcd19a2 = as.Date("2010-07-03"),
    geonexl2 = as.Date("2018-07-03"))
if (!is.null(Wf$test.small.daterange) && Wf$test.small.daterange)
   {Wf$years = year(Wf$date.example)
    Wf$dates = Wf$date.example + (-1:1)}
stopifnot(Wf$date.example %in% Wf$dates)

Wf$y.sat = "Optical_Depth_047"
Wf$features = c(
  # Predictors for modeling.
    "y.sat", "time.sat",
    "y.sat.mean", "y.sat.present",
    "AOD_Uncertainty",
    "cosSZA", "cosVZA", "RelAZ", "Scattering_Angle", "Glint_Angle",
    (if (Wf$satellite.product == "mcd19a2")
        c("Column_WV", "qa_best")))
Wf$window.radius = 5L
  # The windows will be `1 + 2*window.radius` cells on each side.
Wf$pred.round.digits = 5L
  # MCD19A2 AOD has 3 digits of precision, so this is a little more.

# Put the `targets` configuration file and data store in a workflow-
# specific subdirectory of `data.dir`.
workflow.dir = file.path(data.dir, "workflows",
    digest::digest(config$workflow, algo = "murmur32"))
dir.create(workflow.dir, recursive = T, showWarnings = F)
yaml::write_yaml(file = file.path(workflow.dir, "workflow.yaml"),
    config$workflow)
yaml::write_yaml(file = file.path(workflow.dir, "targets.yaml"),
    list(main = list(
        script = "code/targets.R",
        store = file.path(workflow.dir, "targets_store"))))
Sys.setenv(TAR_CONFIG = file.path(workflow.dir, "targets.yaml"))

daily.sat = function(satellite.product = Wf$satellite.product)
  # Whether the satellite product is daily, as opposed to being
  # resolved at a finer unit, such as the second.
    satellite.product != "geonexl2"

multipass.sat = function(satellite.product = Wf$satellite.product)
  # Whether the satellite product has more than one overpass in
  # each file.
    satellite.product == "mcd19a2"

# for reporting purpose
y_var = "diff_AOD"

set.seed(1234)

terra::terraOptions(progress = 0)

# basic AERONET variables to read
vars0 <- c("Date(dd:mm:yyyy)", "Time(hh:mm:ss)", "Day_of_Year","AERONET_Site_Name",
           "Site_Elevation(m)", "Ozone(Dobson)", "NO2(Dobson)",
           "Solar_Zenith_Angle(Degrees)", "Precipitable_Water(cm)")

feature.raster.layers = c(
    "AOD_Uncertainty",
    "cosSZA", "cosVZA", "RelAZ", "Scattering_Angle", "Glint_Angle",
    "Column_WV", "AOD_QA")