# contains functions used in generating the CONUS_AOD paper

agreement.plot.data = \(d)
   {d = melt(d[, .(y.ground, y.ground.pred, y.sat)],
        id = "y.ground", variable.name = "comparison", value.name = "y.pred")
    d[, comparison := factor(comparison, levels = rev(levels(comparison)))]
    d[, in.region := between(y.ground, 0, 1) & between(y.pred, 0, 1)]
    d}

agreement.plot = \(d, base_size = 11)
   {d = agreement.plot.data(d)
    sample.obs = `[`(
        CJ(
            comparison = unique(d$comparison),
            y.ground = seq(0, 1, len = 1000)),
        by = comparison,
        j =
           {ps = with(d[comparison == .BY$comparison],
                empirical.error.envelope(obs = y.ground, y.pred)$par)
            punl(
                y.ground,
                envelope.hi = y.ground + ps[1] + ps[2]*y.ground,
                envelope.lo = y.ground - ps[1] - ps[2]*y.ground)})
    half_line = base_size/2
    ggplot() +
        ggpointdensity::geom_pointdensity(
            data = d[(in.region)], aes(y.ground, y.pred),
            size = 1/4, adjust = 1/2) +
        scale_color_gradient(low = "#dddddd", high = "black") +
        geom_abline(linetype = "dashed", color = "#6666aa") +
        geom_line(data = sample.obs, aes(y.ground, envelope.hi),
            color = "#6666aa") +
        geom_line(data = sample.obs, aes(y.ground, envelope.lo),
            color = "#6666aa") +
        xlab("AERONET AOD") +
        ylab("Satellite-based AOD") +
        guides(color = "none") +
        facet_wrap(vars(comparison), labeller = as_labeller(
          c(`y.sat` = "MAIAC AOD", `y.ground.pred` = "Corrected AOD"))) +
        theme_classic() +
        theme(panel.spacing = unit(1.5, "lines"),
          plot.margin = margin(half_line,
            base_size, half_line, half_line)) +
        coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = F)}

#' Find an "error envelope" with an additive and multiplicative
#' term in the fashion of Dark Target validation
#' https://darktarget.gsfc.nasa.gov/validation
empirical.error.envelope = \(..., target.coverage = 2/3)
   {f = envelope.tester(...)
    optim(c(.05, .1), \(p)
        sum(abs(target.coverage - f(p[1], p[2]))))}

envelope.tester = \(obs, pred, threshold = 0.6)
   {get.coverage = \(slice, add, mult) slice[, mean(
        pred >= obs - add - mult*obs &
        pred <= obs + add + mult*obs)]
    d = data.table(obs, pred)
    d.lo = d[obs <= threshold]
    d.hi = d[obs > threshold]
    \(add, mult) c(
        lo = get.coverage(d.lo, add, mult),
        hi = get.coverage(d.hi, add, mult))}

#' Get observations of PM_{2.5} from the Environmental Protection
#' Agency's (EPA) Air Quality System. The unit is μg/m^3.
#' File source: https://aqs.epa.gov/aqsweb/airdata/download_files.html#Daily
#' Documentation: https://aqs.epa.gov/aqsweb/documents/about_aqs_data.html
get.aqs.obs = function(time.unit, years, grid)
   {assert(time.unit %in% c("hour", "day", "year"))
    aqs.url.root = "https://aqs.epa.gov/aqsweb/airdata"
    parameter.code = 88101L
      # PM_{2.5} from a federally approved reference or equivalent
      # method (FRM/FEM).
    min.days.for.annual = 100L

    assert(file.exists("code/openssl_workaround.conf"))
    Sys.setenv(OPENSSL_CONF = "code/openssl_workaround.conf")
    d = rbindlist(lapply(years, function(the.year)
       {message("AQS ", the.year)
        fname = switch(time.unit,
            hour = sprintf("hourly_%d_%d.zip",
                parameter.code,
                the.year),
            day = sprintf("daily_%d_%d.zip",
                parameter.code,
                the.year),
            year = sprintf("annual_conc_by_monitor_%d.zip",
                the.year))
        d = fread(cmd = paste("unzip -p ", shQuote(download(
            paste0(aqs.url.root, "/", fname),
            file.path("aqs", fname)))))
        setnames(d, str_replace_all(names(d), " ", "."))
        switch(time.unit,
            hour = d[,
                .(
                    time = `-`(
                        lubridate::fast_strptime(
                            paste(Date.GMT, Time.GMT),
                            "%Y-%m-%d %H:%M", lt = F),
                        # Subtract 1 hour, since AQS labels the end rather
                        # than the beginning of each hour spanned.
                        lubridate::hours(1)),
                    lon = Longitude, lat = Latitude,
                    value = Sample.Measurement)],
            day = d[
                Sample.Duration %in% c("24 HOUR", "24-HR BLK AVG") &
                    Event.Type != "Excluded",
                .(
                    time = Date.Local,
                    lon = Longitude, lat = Latitude,
                    value = Arithmetic.Mean)],
            year = d[
                Parameter.Code == parameter.code &
                    Pollutant.Standard == "PM25 Annual 2012" &
                    Required.Day.Count >= min.days.for.annual &
                    Completeness.Indicator == "Y",
                .(
                    time = the.year,
                    lon = Longitude, lat = Latitude,
                    value = as.numeric(Arithmetic.Mean))])}))

    d[, cell := as.integer(terra::cellFromXY(grid,
        convert.crs(cbind(lon, lat), crs.lonlat, terra::crs(grid))))]
    d = d[!is.na(cell)]
    setkey(d, time, cell, lon, lat)
    setcolorder(d)
    d}

get.satellite.vs.aqs = function(time.unit, satellite, aqs)
   {assert(time.unit %in% c("hour", "day", "year"))
    if (time.unit == "year")
        return(get.satellite.vs.aqs.yearly(satellite, aqs))
    message("Merging with AQS")
    d = merge(by = c("time", "cell"),
        satellite[, .(cell, y.sat.old, y.sat.new, time = switch(time.unit,
            day = lubridate::as_date(time.sat, tz = "Etc/GMT+6"),
            hour = lubridate::floor_date(time.sat, "hour")))],
        cbind(aqs, aqs.i = seq_len(nrow(aqs)))[, .(
            aqs.i, cell, y.aqs = value,
            time = switch(time.unit,
                day = time,
                hour = lubridate::floor_date(time, "hour")))])
    # Find the weighted correlation of old and new satellite values
    # with AQS values. A single cell-time can have more than one
    # satellite value, more than one AQS value, or both. We consider
    # all combinations and weight them to sum to 1 per cell-time.
    message("Computing statistics")
    f = \(d)
       {d = copy(d)
        d[, by = .(time, cell), weight := 1 / .N]
        cbind(
            data.table(
                n.aqs.obs = d[, uniqueN(aqs.i)],
                n.cell.times = uniqueN(d[, .(time, cell)])),
            as.data.table(sapply(c("old", "new"), simplify = F, function(sv.type)
                cov.wt(
                    d[, .(get(paste0("y.sat.", sv.type)), y.aqs)],
                    wt = d$weight, cor = T, method = "ML")$cor[1, 2])))}
    if (time.unit == "hour")
        rbind(
            cbind(hour = "all", f(d)),
            d[, keyby = hour(lubridate::with_tz(time, "Etc/GMT+5")),
                f(.SD), .SDcols = colnames(d)])
    else
        f(d)}

get.satellite.vs.aqs.yearly = function(sat, aqs)
   {min.days.in.month = 10L
    min.months = 12L

    # Summarize the satellite predictions into year-cells.
    sat = sat[cell %in% unique(aqs$cell)]
    sat[, date := lubridate::as_date(time.sat, tz = "Etc/GMT+6")]
    sat = sat[, by = .(year(date), month(date), cell), if (uniqueN(date) >= min.days.in.month) .(
        y.sat.old = mean(y.sat.old / 10^Wf$pred.round.digits),
        y.sat.new = mean(y.sat.new / 10^Wf$pred.round.digits))]
    sat = sat[, by = .(time = year, cell), if (uniqueN(month) >= min.months) .(
        y.sat.old = mean(y.sat.old),
        y.sat.new = mean(y.sat.new))]

    aqs = copy(aqs)
    setnames(aqs, "value", "y.aqs")

    merge(by = c("time", "cell"), aqs, sat)[, .(
        n.cell.times = uniqueN(data.table(time, cell)),
        old = cor(y.aqs, y.sat.old),
        new = cor(y.aqs, y.sat.new))]}

get.median.improve.map.data = function(
      y.sat.name, the.satellite, satellite.product, n.workers,
      d, pred.grid, region.shape, satellite.files, model.full)
   {tz = "Etc/GMT+6" # I.e., UTC-6, or Central Standard Time

    # Get the date with the median improvement in MSE.
    date = (d
        [, by = .(date = lubridate::as_date(time.sat, tz = tz)),
            .(improve =
                mean((y.ground - y.ground.pred)^2) /
                mean((y.ground - y.sat)^2))]
        [which.min(improve - median(improve)), date])

    # Get overpasses for this date.
    sat = satellite.files[lubridate::as_date(time, tz = tz) == date]
    sat[, sat.files.ix := .I]
    overpasses = (if (multipass.sat(satellite.product))
        merge(by = "sat.files.ix", sat,
            expand.to.overpasses(sat, sat,
                the.satellite, satellite.product, n.workers)) else
        setnames(cbind(sat, overpass = NA), "time", "time.sat"))

    # Find the local cell indices we'd want to map.
    message("Getting raster cells")
    overpass.local.cells = pblapply(1 : nrow(overpasses), cl = n.workers, \(i)
       {r = with(overpasses[i], read_satellite_raster(
            satellite.product, tile, path, overpass))[[y.sat.name]]
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

    include.quality = "DQF" %in% colnames(d)
    g1 = pred.grid
    raster.vars = c("y.sat.old", "y.sat.new",
        (if (include.quality) "DQF"))
    for (vname in raster.vars)
       {g1[[vname]] = NA_real_
        g1[[vname]][d$cell] = ifelse(
            vname == "DQF" & is.na(d$y.sat.old) & is.na(d$y.sat.new),
            NA_real_,
            d[[vname]])}
    g1 = terra::trim(g1[[raster.vars]])
    g = terra::trim(terra::project(
        g1,
        paste0("epsg:", crs.lonlat),
        res = reproject.res))

    d = melt(
        as.data.table(as.data.frame(g, xy = T, na.rm = F)),
        id.vars = c("x", "y", (if (include.quality) "DQF")))
    setnames(d, c("x", "y"), c("lon", "lat"))
    if (!is.null(limits))
        d = d[
            limits$lon[1] <= lon & lon <= limits$lon[2] &
            limits$lat[1] <= lat & lat <= limits$lat[2]]
    if (!is.null(quantile.cap))
        d[, value := pmin(quantile(value, quantile.cap, na.rm = T), value)]
    if (include.quality)
       d = rbind(
           cbind(d, qualities = "all qualities"),
           cbind(d[DQF != 2], qualities = "excluding worst"))
             # For DQF interpretation, see PDF p. 334 of
             # http://web.archive.org/web/20230915201900if_/https://www.goes-r.gov/products/docs/PUG-L2+-vol5.pdf

    list(
        center.name = bg.sf[as.integer(st_intersects(
            st_as_sf(
                d[, .(
                    mean(c(max(lon), min(lon))),
                    mean(c(max(lat), min(lat))))],
                coords = c(1, 2), crs = crs.lonlat),
            st_transform(bg.sf, crs = crs.lonlat))),]$NAME,
        plot = ggplot() +
            geom_raster(aes(lon, lat, fill = value), data = d) +
            scale_fill_distiller(name = color.scale.name,
                palette = "Spectral", na.value = "transparent") +
            geom_sf(data = bg.sf, fill = NA, size = .1) +
            ggspatial::annotation_scale(data = cbind(
                data.table(variable = "y.sat.old"),
                (if (include.quality) list(qualities = "all qualities")))) +
            facet_grid(
                rows = vars(variable),
                cols = (if (include.quality) vars(qualities)),
                labeller = labeller(variable =
                    c(y.sat.old = "Raw", y.sat.new = "Corrected"))) +
            coord_sf(expand = F,
                xlim = range(d$lon), ylim = range(d$lat)) +
            theme_void() +
            theme(strip.text = element_text(size = 13)))}
