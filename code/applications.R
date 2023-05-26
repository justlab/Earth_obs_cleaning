preds.ctdeep = \(output.dir, pred.grid, region.shape)
  # Save predictions for Michael Geigert at the Connecticut
  # Department of Energy & Environmental Protection (DEEP).
  # This takes about an hour on Coco.
   {dates = seq(as.Date("2021-07-24"), as.Date("2021-07-30"), by = 1)
    tz = "America/New_York"
    corner = data.table(lon = -90, lat = 36)

    message("Getting cell indices")
    corner = convert.crs(corner, crs.lonlat, terra::crs(pred.grid))
    gridrows = (as.data.table(as.data.frame(pred.grid, cells = T, xy = T))
       [x >= corner$x & y >= corner$y])
    setkey(gridrows, cell)

    for (date.i in seq_along(dates))
       {message("~~~~~ ", dates[date.i])
        dt.start = lubridate::as_datetime(dates[date.i], tz = tz)
        dt.end = lubridate::as_datetime(dates[date.i] + 1, tz = tz) - 1
        p = new.preds(dt.start, dt.end, gridrows$cell)[,
            .(time.sat, cell, y.sat.old, y.sat.new)]
        # Reduce the predictions to the most common time per tile.
        p[, tile := gridrows[.(p$cell), tile]]
        p = p[, by = tile,
           {ts = unique(time.sat)
            .SD[time.sat == ts[which.max(tabulate(match(time.sat, unique(ts))))]]}]
        # Make a raster and save it.
        g = pred.grid
        g$aod_original = NA_real_
        g$aod_original[p$cell] = p$y.sat.old
        g$aod_corrected = NA_real_
        g$aod_corrected[p$cell] = p$y.sat.new
        terra::writeRaster(
            terra::trim(g[[c("aod_original", "aod_corrected")]]),
            file.path(output.dir, paste0(dates[date.i], ".tiff")),
            filetype = "GTiff", gdal = c("COMPRESS=DEFLATE"),
            overwrite = T)}}
