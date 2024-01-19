Sys.setenv(RENV_CONFIG_STARTUP_QUIET = "TRUE")
renv::load()

suppressPackageStartupMessages(
   {library(data.table)
    library(targets)})

data.dir = "/data"
stopifnot(dir.exists(data.dir))
aodc.dir = file.path(data.dir, "aodc")
writing.out.dir = file.path(data.dir, "writing")
dir.create(writing.out.dir, showWarnings = F)
intermediate.path = function(...)
   file.path(data.dir, 'intermediate', ...)
dir.create(intermediate.path(), showWarnings = F)
download = function(from, to, ...)
    download.update.meta(from, file.path(data.dir, "downloads"), to, ...)
satellite_hdf_root = file.path(data.dir, 'earthdata')
dir.create(satellite_hdf_root, showWarnings = F)

config = yaml::read_yaml(file.path(data.dir, "config.yaml"))

n.workers = config$n.workers
if (is.null(n.workers))
    n.workers = max(1, parallel::detectCores() - 2)

# The most general workflow-defining variables.
Wf = as.environment(config$workflow)
stopifnot(
    Wf$outcome == "aod" &&
        Wf$satellite.product == "mcd19a2" &&
        Wf$satellite %in% c("terra", "aqua") &&
        Wf$ground.product == "aeronet" ||
    Wf$outcome == "aod" &&
        Wf$satellite.product == "aodc" &&
        Wf$satellite == "goes16" &&
        Wf$ground.product == "aeronet" ||
    Wf$outcome == "no2" &&
        Wf$satellite.product == "tropomi" &&
        Wf$satellite == "sentinel5p" &&
        Wf$ground.product == "pandonia")
stopifnot(is.character(Wf$region))

# True configurability of other `Wf` items is not implemented now
# but might be added later.

Wf$dates = switch(Wf$satellite.product,
    mcd19a2 = c("2000-01-01", "2022-12-31"),
    aodc = c("2020-10-25", "2023-12-31"),
      # "The land spectral surface relationship was updated on October 24, 2020."
      # We start a day later in case of off-by-one issues.
      # http://web.archive.org/web/20240112172248/https://www.star.nesdis.noaa.gov/atmospheric-composition-training/documents/GOES-16_ABI_L2_AOD_Provisional_ReadMe_v3.pdf
    tropomi = c("2021-07-01", "2022-12-31"))
      # The former is the first day of the latest version of the
      # TROPOMI nitrogen-dioxide product.
Wf$dates = seq(as.Date(Wf$dates[1]), as.Date(Wf$dates[2]), by = 1)
Wf$time.example = switch(Wf$satellite.product,
  # This should be a date (or datetime, for more temporally resolved
  # products) for which the satellite data of interest exists on all
  # tiles.
    mcd19a2 = as.Date("2010-07-03"),
    aodc = lubridate::as_datetime("2020-11-21 20:27:30 UTC"),
    tropomi = as.Date("2021-08-01"))
if (!is.null(Wf$test.small.daterange) && Wf$test.small.daterange)
    Wf$dates = lubridate::as.date(Wf$time.example) + (-1:1)
Wf$years = sort(unique(year(Wf$dates)))

Wf$y.sat = switch(Wf$satellite.product,
    mcd19a2 = "Optical_Depth_047",
    aodc = "AOD")
Wf$features = c(
  # Predictors for modeling.
    "y.sat", "time.sat",
    "y.sat.mean", "y.sat.present",
    switch(Wf$satellite.product,
        mcd19a2 = c(
            "qa_best", "AOD_Uncertainty",
            "cosSZA", "cosVZA", "RelAZ", "Scattering_Angle", "Glint_Angle",
            "Column_WV"),
        aodc = c(
            "DQF", "seconds.since.midnight")))
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
    satellite.product != "aodc"

multipass.sat = function(satellite.product = Wf$satellite.product)
  # Whether the satellite product has more than one overpass in
  # each file.
    satellite.product == "mcd19a2"

# for reporting purpose
y_var = "diff_AOD"

terra::terraOptions(progress = 0)

# basic AERONET variables to read
vars0 <- c("Date(dd:mm:yyyy)", "Time(hh:mm:ss)", "Day_of_Year","AERONET_Site_Name",
           "Site_Elevation(m)", "Ozone(Dobson)", "NO2(Dobson)",
           "Solar_Zenith_Angle(Degrees)", "Precipitable_Water(cm)")

feature.raster.layers = switch(Wf$satellite.product,
    mcd19a2 = c(
        "AOD_Uncertainty",
        "cosSZA", "cosVZA", "RelAZ", "Scattering_Angle", "Glint_Angle",
        "Column_WV", "AOD_QA"),
    aodc = c(
        "DQF"))
