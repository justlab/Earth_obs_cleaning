# contains functions used in generating the CONUS_AOD paper

#' Print a nice html table using the `DT` package
#'
dtwrapper <- function(dt){
  DT::datatable(dt,
                rownames = FALSE,
                fillContainer = FALSE,
                autoHideNavigation = TRUE,
                options = list(pageLength = 50,
                               autoWidth = TRUE))
}

#' Performance metrics for MCD19A2 AOD
#'
performance_aod <- function(dt, ...){
  performance_metrics(dt,
                      ground_truth = "AOD_470nm",
                      eo_raw = "MCD19_AOD_470nm",
                      eo_pred = "aodhat")
}

#' Get observations of PM_{2.5} from the Environmental Protection
#' Agency's (EPA) Air Quality System. The unit is μg/m^3.
#' File source: https://aqs.epa.gov/aqsweb/airdata/download_files.html#Daily
#' Documentation: https://aqs.epa.gov/aqsweb/documents/about_aqs_data.html
get_ground_obs = function(years, grid)
   {aqs.url.root = "https://aqs.epa.gov/aqsweb/airdata"
    parameter.code = 88101L
      # PM_{2.5} from a federally approved reference or equivalent
      # method (FRM/FEM).

    assert(file.exists("openssl_workaround.conf"))
    Sys.setenv(OPENSSL_CONF = "openssl_workaround.conf")
    d = rbindlist(lapply(years, function(the.year)
       {fname = sprintf("daily_%d_%d.zip", parameter.code, the.year)
        d = fread(
            cmd = paste("unzip -p ", shQuote(download(
                paste0(aqs.url.root, "/", fname),
                file.path("aqs", fname)))),
            select = c(
               "Date Local", "Longitude", "Latitude",
               "Event Type", "Sample Duration", "Arithmetic Mean"))
        setnames(d, str_replace_all(names(d), " ", "."))
        d[
            Sample.Duration %in% c("24 HOUR", "24-HR BLK AVG") &
                Event.Type != "Excluded",
            .(date = Date.Local, lon = Longitude, lat = Latitude,
               value = Arithmetic.Mean)]}))

    d[, cell := as.integer(terra::cellFromXY(grid,
        convert.crs(cbind(lon, lat), crs.lonlat, terra::crs(grid))))]
    d = d[!is.na(cell)]
    setkey(d, date, cell, lon, lat)
    setcolorder(d)
    d}

#' Get all satellite observations at cells that ever contain AQS
#' monitors.
satellite_at_aqs_sites = function(region, years, sat, ground)
    rbindlist(lapply(
        as.character(outer(1:12, years, sprintf,
            fmt = "pred_out_%d_%d_%s_%s", sat, region)),
        function(tname){
          satellite = tar_read_raw(tname)
          if (nrow(satellite) == 0) NULL else satellite[
            cell %in% unique(ground$cell),
            .(
              date = pred_date,
              overpass, cell,
              satellite_value_old = value_old,
              satellite_value_new = value_new)]
        }
    ))

satellite_vs_ground = function(satellite, ground)
   {d = merge(
       satellite,
       ground[, .(date, cell, ground_value = value)],
       by = c("date", "cell"))
    # Find the weighted correlation of old and new satellite values
    # with ground values.
    d[, by = .(date, cell), weight := 1 / .N]
    out = lapply(c("old", "new"), function(sv_type)
        cov.wt(
            d[, .(get(paste0("satellite_value_", sv_type)), ground_value)],
            wt = d$weight, cor = T, method = "ML")$cor[1, 2])
    names(out) <- c("old", "new")
    out[["nobs"]] = d[, .N]
    out[["celldays"]] = d[, sum(weight)]
    out
    }
