new.preds = function(dt.start, dt.end, cells = NULL, targets = NULL)
  # Make a data table of new predictions for every non-missing
  # satellite value that occurs in the given spacetime chunk.
  # - `dt.start` and `dt.end` specify a range of datetimes to make
  #   predictions for, as `POSIXct` objects.
  # - `cells`, if provided, should be a vector of cell numbers of
  #   `pred_grid`. Otherwise, we use all cells.
   {if (is.null(targets))
       {message("Reading targets")
        grid = tar_read(pred_grid)
        sat = tar_read(satellite_hdf_files)
        model = tar_read(full_model)}
    else
       {grid = targets[[1]]
        sat = targets[[2]]
        model = targets[[3]]}

    message("Checking inputs")
    for (dt.v in list(dt.start, dt.end))
        assert("POSIXct" %in% class(dt.v))
    dt.start = lubridate::with_tz(dt.start, "UTC")
    dt.end = lubridate::with_tz(dt.end, "UTC")
    if (is.null(cells))
        cells = which(!is.na(drop(grid$tile[])))
    if (is.double(cells))
        cells = as.integer(cells)
    assert(is.integer(cells) && all(
        1L <= cells & cells <= terra::ncell(grid)))
    message("Cells: ", scales::comma(length(cells)))

    message("Selecting satellite files")
    # Keep only enough `sat` files to cover the requested time range.
    sat = sat[
        if (class(time) == "Date")
            lubridate::as_date(dt.start) - 1 <= time &
                time <= lubridate::as_date(dt.end) + 1
        else
            dt.start <= time & time <= dt.end]
    # And keep only files for the necessary tiles.
    sat = sat[tile %in% unique(grid$tile[cells][[1]])]
    sat[, sat.files.ix := .I]

    d = expand.to.overpasses(
        sat[, .(time.sat = time, sat.files.ix)],
        sat,
        Wf$satellite, Wf$satellite.product, config$n.workers)
    if (!nrow(d))
        return(data.table())

    message("Expanding to one row per cell-time")
    d = d[dt.start <= time.sat & time.sat <= dt.end]
    d[, tile := sat[d$sat.files.ix, tile]]
    d = merge(d, cbind(as.data.table(grid[cells]), cell = cells),
        by = "tile", allow.cartesian = T)
    message(sprintf("Result: %s rows", scales::comma(nrow(d))))

    d = get.predictors(d, sat,
        Wf$satellite.product, Wf$y.sat,
        Wf$features, Wf$window.radius,
        config$n.workers)
    message(sprintf("Result: %s rows", scales::comma(nrow(d))))
    if (!nrow(d))
        return(data.table())

    message("Making new predictions")
    simplify.dt.for.xgboost(d)
    d[, y.sat.new := y.sat - predict(
        xgb.load.raw(model$model),
        as.matrix(d[, mget(Wf$features)]))]
    d[, time.sat := lubridate::as_datetime(time.sat)]
    setnames(d, "y.sat", "y.sat.old")

    message("Sorting")
    setkey(d, time.sat, cell)
    d}

new.preds.compact = function(...)
   {r = function(v) as.integer(round(10^Wf$pred.round.digits * v))
    d = new.preds(...)
    if (!nrow(d))
        return(data.table())
    d[, .(
        time.sat,
        cell,
        y.sat.old = r(y.sat.old),
        y.sat.new = r(y.sat.new))]}
