library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Globals ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Region
aoiname = "conus"

satellite_aod_tiles = list(
    # Selected by hand in QGIS with this shapefile:
    # http://web.archive.org/web/2022/http://book.ecosens.org/wp-content/uploads/2016/06/modis_grid.zip
    conus = rbind(
        data.table(h = 8:13, v = 4),
        data.table(h = 8:12, v = 5),
        data.table(h = 8:10, v = 6))[, sprintf("h%02dv%02d", h, v)])

# which MODIS platforms to use
sats = c("terra", "aqua")

# out_dir = paste0("/data-coco/mcd19/fst/conus/", this_year, "_AOD_newVar/")
# if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# for reporting purpose
y_var = "diff_AOD"

set.seed(1234)

# basic AERONET variables to read
vars0 <- c("Date(dd:mm:yyyy)", "Time(hh:mm:ss)", "Day_of_Year","AERONET_Site_Name",
           "Site_Elevation(m)", "Ozone(Dobson)", "NO2(Dobson)",
           "Solar_Zenith_Angle(Degrees)", "Precipitable_Water(cm)")

crs_sinu = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
# # temporary AOI for testing prediction table preparation
# library(sf)
# dc_wgs = st_bbox(c(xmin = -77.78, xmax = -76.45, ymax = 39.27, ymin = 38.52), crs = 4326)
# dc_sinu = st_bbox(st_transform(st_as_sfc(dc_wgs), crs_sinu))
# rm(dc_wgs)
