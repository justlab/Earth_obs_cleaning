# contains functions used in generating the CONUS_AOD paper

#' Try out empirical error envelope parameters
#'
empirical.error.envelope <- function(dt,
                                     obs = "y.ground",
                                     pred = "y.sat",
                                     add.term = seq(0.01, 0.05, by = 0.0025),
                                     mult.term = seq(0.125, 0.2, by = 0.0025),
                                     threshold = 0.6,
                                     coverage = 2/3){
  envelope.terms <- CJ(add.term, mult.term)
  coverage.percent <- function(i, dt) {
    temp <- unlist(i)
    with(dt, mean(get(pred) > get(obs) - temp[["add.term"]] - temp[["mult.term"]]*get(obs) &
                            get(pred) < get(obs) + temp[["add.term"]] + temp[["mult.term"]]*get(obs)))
  }
  envelope.terms[, cov.low := coverage.percent(.SD, dt[get(obs) <= threshold]), by = row.names(envelope.terms)]
  envelope.terms[, cov.high := coverage.percent(.SD, dt[get(obs) > threshold]), by = row.names(envelope.terms)]
  envelope.terms[, cov.overall := coverage.percent(.SD, dt), by = row.names(envelope.terms)]
  # return the parameter combo with the smallest abs diff from desired coverage
  envelope.terms <- envelope.terms[order(abs(coverage - cov.low) + abs(coverage - cov.high), decreasing = FALSE),]
  envelope.terms
}

#' Get observations of PM_{2.5} from the Environmental Protection
#' Agency's (EPA) Air Quality System. The unit is μg/m^3.
#' File source: https://aqs.epa.gov/aqsweb/airdata/download_files.html#Daily
#' Documentation: https://aqs.epa.gov/aqsweb/documents/about_aqs_data.html
get.aqs.obs = function(years, grid)
   {aqs.url.root = "https://aqs.epa.gov/aqsweb/airdata"
    parameter.code = 88101L
      # PM_{2.5} from a federally approved reference or equivalent
      # method (FRM/FEM).

    assert(file.exists("code/openssl_workaround.conf"))
    Sys.setenv(OPENSSL_CONF = "code/openssl_workaround.conf")
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

get.satellite.vs.aqs = function(satellite, aqs)
   {message("Merging with AQS")
    d = merge(
       satellite[, .(cell, y.sat.old, y.sat.new,
           date = lubridate::as_date(time.sat, tz = "Etc/GMT+6"))],
       aqs[, .(date, cell, y.aqs = value)],
       by = c("date", "cell"))
    # Find the weighted correlation of old and new satellite values
    # with AQS values. A single cell-day can have more than one
    # satellite value, more than one AQS value, or both. We consider
    # all combinations and weight them to sum to 1 per cell-day.
    message("Computing statistics")
    d[, by = .(date, cell), weight := 1 / .N]
    out = lapply(c("old", "new"), function(sv_type)
        cov.wt(
            d[, .(get(paste0("y.sat.", sv_type)), y.aqs)],
            wt = d$weight, cor = T, method = "ML")$cor[1, 2])
    names(out) <- c("old", "new")
    out[["nobs"]] = d[, .N]
    out[["celldays"]] = d[, sum(weight)]
    out
    }

get.median.mse.map.data = function(
      y.sat.name, the.satellite, satellite.product, n.workers,
      d, pred.grid, region.shape, satellite.files, model.full)
   {# Get the date with the median MSE.
    date = (d
        [, by = .(date = lubridate::as_date(time.sat, tz = "UTC")),
            .(mse = mean((y.ground - y.ground.pred)^2))]
        [which.min(mse - median(mse)), date])

    # Get overpasses for this data.
    sat = satellite.files[time == date]
    sat[, sat.files.ix := .I]
    overpasses = merge(by = "sat.files.ix", sat,
        expand.to.overpasses(sat, sat,
            the.satellite, satellite.product, n.workers))

    # Find the local cell indices we'd want to map.
    message("Getting raster cells")
    overpass.local.cells = pblapply(1 : nrow(overpasses), cl = n.workers, \(i)
       {r = with(overpasses[i], read_satellite_raster(
            Wf$satellite.product, tile, path, overpass))[[y.sat.name]]
        # Get only cells with non-missing values.
        d = as.data.table(as.data.frame(r,
            na.rm = T, cell = T, xy = T))
        if (!nrow(d))
            return(integer())
        # Get only cells in the study area.
        d[in.sf(x, y, terra::crs(r), region.shape), cell]})

    # Choose the tile-time with the most non-missing values.
    oi = which.max(sapply(overpass.local.cells, length))
    # Get the prediction-grid cell for each local cell.
    message("Getting cell indices")
    cells = (as.data.table(as.data.frame(pred.grid, cells = T))
        [tile == overpasses[oi, tile]]
        [cell.local %in% overpass.local.cells[[oi]], cell])
    assert(length(cells) == length(overpass.local.cells[[oi]]))

    # Make the predictions and return.
    list(
        date = date,
        time.sat = overpasses[oi, time.sat],
        tile = overpasses[oi, tile],
        pred = new.preds(
            overpasses[oi, time.sat],
            overpasses[oi, time.sat],
            cells,
            targets = list(pred.grid, satellite.files, model.full)))}

get.baltimore.map.data = function(pred.grid, region.shape, satellite.files, model.full)
   {tiles = c("h11v05", "h12v05")
    dt = lubridate::as_datetime("2015-06-10T15:40:00Z")

    message("Getting cell indices")
    cells = (as.data.table(as.data.frame(pred.grid, cells = T, xy = T))
       [tile %in% tiles]
       [in.sf(x, y, terra::crs(pred.grid), region.shape), cell])

    new.preds(
        dt, dt, cells,
        targets = list(pred.grid, satellite.files, model.full))}

pred.map = function(
        d, pred.grid, bg.sf, color.scale.name,
        limits = NULL, quantile.cap = NULL)
   {reproject.res = .009

    g1 = pred.grid
    g1$y.sat.old = NA_real_
    g1$y.sat.new = NA_real_
    g1$y.sat.old[d$cell] = d$y.sat.old
    g1$y.sat.new[d$cell] = d$y.sat.new
    g1 = terra::trim(g1[[c("y.sat.old", "y.sat.new")]])
    g = terra::trim(terra::project(
        g1,
        paste0("epsg:", crs.lonlat),
        res = reproject.res))
    #print(data.table(
    #    Grid = c("Original", "Reprojected"),
    #    "Non-NA cells" = scales::comma(c(sum(!is.na(g1[])), sum(!is.na(g[]))))))

    d = melt(
        as.data.table(as.data.frame(g, xy = T, na.rm = F)),
        id.vars = c("x", "y"))
    setnames(d, c("x", "y"), c("lon", "lat"))
    if (!is.null(limits))
        d = d[
            limits$lon[1] <= lon & lon <= limits$lon[2] &
            limits$lat[1] <= lat & lat <= limits$lat[2]]
    if (!is.null(quantile.cap))
        d[, value := pmin(quantile(value, quantile.cap, na.rm = T), value)]

    ggplot() +
        geom_raster(aes(lon, lat, fill = value), data = d) +
        scale_fill_distiller(name = color.scale.name,
            palette = "Spectral", na.value = "transparent") +
        geom_sf(data = bg.sf, fill = NA, size = .1) +
        facet_grid(rows = "variable", labeller = labeller(variable =
            c(y.sat.old = "Original", y.sat.new = "Corrected"))) +
        coord_sf(expand = F,
            xlim = range(d$lon), ylim = range(d$lat)) +
        theme_void()}
