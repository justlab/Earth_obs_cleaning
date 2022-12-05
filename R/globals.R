library(data.table)

data.dir = Sys.getenv("EARTH_OBS_CLEANING_DATA_DIR")
stopifnot(dir.exists(data.dir))
geonexl2.dir = file.path(data.dir, "geonexl2")
intermediate.path = function(...)
   file.path(data.dir, 'intermediate', ...)
dir.create(intermediate.path(), showWarnings = F)
download = function(from, to, ...)
    download.update.meta(from, file.path(data.dir, "downloads"), to, ...)
satellite_hdf_root = file.path(data.dir, 'earthdata')
dir.create(satellite_hdf_root, showWarnings = F)

n.workers = Sys.getenv("EARTH_OBS_CLEANING_NTHREADS")
n.workers = (if (n.workers == "")
    max(1, parallel::detectCores() - 2) else
    as.integer(n.workers))

# The most general workflow-defining variables.
Wf = stringr::str_match_all(
   Sys.getenv("EARTH_OBS_CLEANING_WORKFLOW"),
   "([^ =]+)=([^ =]+)")[[1]]
Wf = as.environment(`names<-`(as.list(Wf[,3]), Wf[,2]))
stopifnot(Wf$outcome %in% c("aod"))
stopifnot(
    Wf$satellite.product == "mcd19a2" &&
       Wf$satellite %in% c("terra", "aqua") ||
    Wf$satellite.product == "geonexl2" &&
       Wf$satellite == "goes16")
stopifnot(Wf$ground.product %in% c("aeronet"))
stopifnot(Wf$region %in% c("conus"))

daily.sat = function(satellite.product = Wf$satellite.product)
  # Whether the satellite product is daily, as opposed to being
  # resolved at a finer unit, such as the second.
    satellite.product != "geonexl2"

multipass.sat = function(satellite.product = Wf$satellite.product)
  # Whether the satellite product has more than one overpass in
  # each file.
    satellite.product == "mcd19a2"

satellite_aod_tiles = list(
    # Selected by hand in QGIS with this shapefile:
    # http://web.archive.org/web/2022/http://book.ecosens.org/wp-content/uploads/2016/06/modis_grid.zip
    conus = rbind(
        data.table(h = 8:13, v = 4),
        data.table(h = 8:12, v = 5),
        data.table(h = 8:10, v = 6))[, sprintf("h%02dv%02d", h, v)])

# for reporting purpose
y_var = "diff_AOD"

set.seed(1234)

terra::terraOptions(progress = 0)

# basic AERONET variables to read
vars0 <- c("Date(dd:mm:yyyy)", "Time(hh:mm:ss)", "Day_of_Year","AERONET_Site_Name",
           "Site_Elevation(m)", "Ozone(Dobson)", "NO2(Dobson)",
           "Solar_Zenith_Angle(Degrees)", "Precipitable_Water(cm)")

crs_sinu = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
