## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## * Common modeling code
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

simplify.dt.for.xgboost = function(d)
    if ("time.sat" %in% colnames(d))
        d[, time.sat := as.double(time.sat)]

#' internal function to get threads to use
#'
get.threads <- function(){
  parallel::detectCores()%/%2
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## * CV
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Do initial CV with hyperparameter selection with DART and rank features.
#'
#' Must provide either stn_var or day_var for binning folds.
#'
#' @param data data.frame of observations with DV and IVs
#' @param y_var column name to predict
#' @param features character vector of column names to train on
#' @param stn_var if binning by stations, provide the station column name
#' @param day_var if binning by days, provide the day column name
# Depends on functions in xgboost_cv_RFE.R
cv_dart <- function(
  data,
  y_var,
  features,
  k_fold = 5,
  n_rounds = 100,
  stn_var = NULL,
  day_var = NULL,
  progress = TRUE
){
  xgb_threads <- get.threads()

  mDT <- data.table::copy(setDT(data)) # needed for drake, uncertain about targets
  simplify.dt.for.xgboost(mDT)
  if(!is.null(day_var)) {
    mDT[, dayint := as.integer(get(day_var))]
    by_var = "day"
  }
  if(!is.null(stn_var)) {
    mDT[, stn := get(stn_var)]
    by_var = "stn"
  }
  if (!is.null(day_var) & !is.null(stn_var))
    stop("Please provide either day_var or stn_var.")

  # some internal variables names
  y_formula <- as.formula(paste(y_var, "~."))
  y_var_pred <- paste0(y_var, "_pred") # name of the predicted y
  y_var_pred_whole <- paste0(y_var, "_pred_whole") # name of the predicted y

  # bin data into specified number of folds
  temp <- prepare.bin(mDT, by_var = by_var, k_fold = k_fold)
  mDT <- temp$data
  bin_list <- temp$bin # list of values of selected variable (stn or day) by fold
  rm(temp)
  index.fs.list <-  list(
    # return the observations in each fold
    stn = function(x, bin_list) mDT[stn%in%bin_list[[x]], which = TRUE],
    day = function(x, bin_list) mDT[dayint%in%bin_list[[x]], which = TRUE]
  )

  index_train <- index_test <- list()
  for (k in 1:k_fold){
    index_test[[k]]  <- index.fs.list[[by_var]](k, bin_list)
    index_train[[k]] <- -index_test[[k]]
  }

  # run k-fold cv and record SHAP matrix, predicted value, overall rmse...
  if (!y_var%in%features) features <- c(features, y_var)

  message("Run k-fold cv \n")
  cv_results <- run.k.fold.cv(k_fold = k_fold,
                              dataXY_df = mDT[, ..features],
                              by_var = by_var,
                              n_rounds = n_rounds,
                              y_var = y_var,
                              progress = progress,
                              seed = 1234,
                              index_train = index_train,
                              index_test = index_test,
                              xgb_threads = xgb_threads)

  mDT <- cbind(mDT, cv_results$y_pred_dt)
  mDT[, time.sat := lubridate::as_datetime(time.sat)]
  mDT[, y.ground.pred := y.sat - y.diff_pred]

  # output
  c(list(by_var = by_var, bin_list = bin_list, mDT_wPred = mDT), cv_results)
}

#' core cv function what uses \code{xgboost.dart.cvtune} to run cv for one fold.
#'
#' internal function.
#' @inheritParams run.k.fold.cv.rfe.wrap
#' @param dataXY_df dataset with only X and Y
#' @param index_train index for training
#' @param index_test index for testing
#' @param xgb_threads xgb_threads passed from parent function
#' @param by_var string of "stn" or "day" passed from parent function
#' @param progress whether to show progress bars, default TRUE
#' @param ... other arguments
#'
run.k.fold.cv <- function(k_fold, dataXY_df, y_var,
                          index_train, index_test, xgb_threads, by_var, n_rounds,
                          progress = TRUE, seed = 1234, ...){
  y_var_pred <- paste0(y_var, "_pred") # name of the predicted y
  Y <-  dataXY_df[, ..y_var]
  data_X <- dataXY_df[, -..y_var]
  # prepare dataset
  n_row <- nrow(data_X)
  n_col <- ncol(data_X)


  # output dataset
  y_pred_dt <- data.table(y_pred = rep(NA_real_, n_row))
  shap_score <- as.data.table(matrix(rep(NA_real_, n_row*n_col), ncol = n_col))
  names(shap_score) <- names(data_X)
  xgb_param_list1 <- xgb_param_list2 <- list()
  BIAS0 <- rep(NA_real_, k_fold)
  # loop for each fold

  if(!(is.null(seed) || is.na(seed))) set.seed(seed)
  for (i in 1:k_fold){
    cat('number of obs in testing fold', i, 'is:', nrow(data_X[index_test[[i]], ]), '\n')

    rsxgb0 <- xgboost.dart.cvtune(
      # by default, gives 100 rounds, and it is enough by experience
      n.rounds = n_rounds,
      d = dataXY_df[index_train[[i]],], dv = y_var, ivs = colnames(data_X),
      objective = "reg:squarederror",
      eval_metric = "rmse",
      progress = progress, nthread = xgb_threads)
    # manually select and store some params
    xgb_param_dart <- c(rsxgb0$model$params[c(1,2,4, 6:11)], nrounds = rsxgb0$model$niter)
    xgbmod <- rsxgb0$model

    # fit model
    y_pred_dt[index_test[[i]], y_pred:= rsxgb0$pred.fun(dataXY_df[index_test[[i]],])] # record the prediction
    # predicted SHAP
    shap_pred <- as.data.table(rsxgb0$pred.fun(dataXY_df[index_test[[i]],],
                                               predcontrib = TRUE, approxcontrib = FALSE))

    BIAS0[i] <- first(shap_pred$BIAS)
    shap_pred[, BIAS := NULL]
    shap_score[index_test[[i]],] <- shap_pred # record the SHAP values (for whole model)

    xgb_param_list1[[paste0(by_var, i)]] <- unlist(xgb_param_dart)
    xgb_param_list2[[paste0(by_var, i)]] <- xgb_param_dart
  }
  cat("loop finished\n")
  rmse_all_folds <- sqrt(mean((y_pred_dt$y_pred - dataXY_df[[y_var]])^2))
  # features ranked by SHAP
  features_rank_full_model <- names(colMeans(abs(shap_score))[order(colMeans(abs(shap_score)),
                                                                    decreasing = TRUE)])
  setnames(y_pred_dt, "y_pred", y_var_pred)
  return(list(rmse_all_folds = rmse_all_folds,
              features_rank_full_model = features_rank_full_model,
              shap_score = shap_score,
              y_pred_dt = y_pred_dt,
              BIAS = BIAS0,
              xgb_param_list1 = xgb_param_list1,
              xgb_param_list2 = xgb_param_list2))
}

#' bin the data into the specified number of folds by day or by station.
#'
#' Internal function for cv
#'
#' @param mDT dataset
#' @param by_var must be either stn or day
#' @param k_fold k fold
#' @return list of a dataset with group assignment and a `bin_list` listing members of each bin
#'
prepare.bin <- function(mDT, by_var, k_fold){
  if(!by_var %in% c('stn', 'day')) stop('choose between by_var = "stn" or "day"')
  ss <- unique(mDT[, .(group_count = .N), by = by_var])
  ss <- ss[order(group_count, decreasing = TRUE),]

  ## use helper.pack.bins is better
  ss$bin_stn <- helper.pack.bins(ss[, group_count], k_fold)
  bin_list <- split(ss[[by_var]], ss$bin_stn)

  setkeyv(mDT, by_var)
  setkeyv(ss, by_var)
  mDT <- mDT[ss, on = by_var]
  mDT[, bin_stn := as.factor(bin_stn)]

  message('Obs in each stn_bin and total obs:')
  print(ss[, .(count = sum(group_count)), by = bin_stn][order(bin_stn)])

  return(list(data = mDT, bin = bin_list))
}

#' assign groups numbers to each station by their number of observations.
#'
#' Internal function for cv
#'
#' @param object.sizes a vector of object size
#' @param n.bins number of bins (e.g. 5 folders )
#' @return a vector of group assignment
#' @examples helper.pack.bins(seq(1:10), n.bins = 4)
helper.pack.bins = function(object.sizes, n.bins){
  # by Kodi
  # Greedily pack objects with sizes given by `object.sizes`
  # into `n.bins` bins. Return a vector of the bin indices.
  d = data.table::data.table(size = as.integer(object.sizes),
                             bin = NA_integer_)
  d[, oi := .I]
  # Shuffle the input so ties on object size are broken randomly.
  d = d[sample.int(nrow(d))]
  bin.sizes = rep(0L, n.bins)
  while (anyNA(d$bin)){
    the.oi = d[is.na(bin)][which.max(size), oi]
    bsi = which.min(bin.sizes)
    d[oi == the.oi, bin := bsi]
    bin.sizes[bsi] = bin.sizes[bsi] + d[oi == the.oi, size]
  }
  d[order(oi), bin]
}

#' Summarize CV statistics.
#'
#' @param d The data table `mDT_wPred` returned by `cv_dart`.
get.cv.summary = function(d)
    rbind(
        d[, c(list(Year = "all"), eval(performance.j))],
        d[, keyby = .(Year = year(time.sat)), eval(performance.j)])

performance.j = quote(
   {mse = function(x, y) mean((x - y)^2)
    sdn = function(x) sqrt(mse(x, mean(x)))
    tqs = unname(quantile(
        abs(as.numeric(difftime(units = "mins",
            time.sat, time.ground))),
        c(.25, .5, .75)))
    list(
        "Cases" = .N,
        "Sites" = length(unique(site)),
        "RMSE, raw" = sqrt(mse(y.ground, y.sat)),
        "RMSE, corrected" = sqrt(mse(y.ground, y.ground.pred)),
        "Proportion of raw MSE" = mse(y.ground, y.ground.pred) / mse(y.ground, y.sat),
        "Median, ground" = median(y.ground),
        "SD, ground" = sdn(y.ground),
        "SD, raw" = sdn(y.sat),
        "SD, corrected" = sdn(y.ground.pred),
        "Bias, raw" = mean(y.sat - y.ground),
        "Bias, corrected" = mean(y.ground.pred - y.ground),
        "Correlation, raw" = cor(y.sat, y.ground),
        "Correlation, corrected" = cor(y.ground.pred, y.ground),
        "Time difference, Q1" = tqs[1],
        "Time difference, median" = tqs[2],
        "Time difference, Q3" = tqs[3])})

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## * New predictions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Train a full model using DART and 2-fold CV to select hyperparameters from
#' maximin Latin hypercube sample
#'
dart_full <- function(
  data_train,
  y_var,
  features,
  n_rounds = 100,
  progress = TRUE
){
  xgb_threads <- get.threads()
  data_train = data_train[, c(features, y_var), with = F]
  simplify.dt.for.xgboost(data_train)
  xdc_out <- xgboost.dart.cvtune(
    # by default, gives 100 rounds, and it is enough by experience
    n.rounds = n_rounds,
    d = data_train,
    dv = y_var,
    ivs = features,
    progress = progress,
    nthread = xgb_threads)
  # record the prediction
  preds = xdc_out$pred.fun(data_train)
  #rmse_full <- sqrt(mean((preds - data_train[[y_var]])^2))

  # SHAP
  shap_pred <- as.data.table(xdc_out$pred.fun(data_train, predcontrib = TRUE, approxcontrib = FALSE))
  shap_bias <- first(shap_pred$BIAS)
  shap_pred[, BIAS := NULL]

  # features ranked by SHAP
  mean_shaps = colMeans(abs(shap_pred))
  ## feature names only:
  # features_rank_full_model <- names(mean_shaps)[order(mean_shaps, decreasing = TRUE)]
  # named vector:
  features_rank_full_model = mean_shaps[order(mean_shaps, decreasing = T)]

  list(features_rank_full_model = features_rank_full_model,
       y_preds = preds,
       shap_pred = shap_pred,
       shap_bias = shap_bias,
       model = xgb.save.raw(raw_format = "ubj", xdc_out$model))}

new.preds = function(dt.start, dt.end, cells = NULL, targets = NULL)
  # Make a data table of new predictions for every non-missing
  # satellite value that occurs in the given spacetime chunk.
  # - `dt.start` and `dt.end` specify a range of datetimes to make
  #   predictions for, as `POSIXct` objects.
  # - `cells`, if provided, should be a vector of cell numbers of
  #   `pred.grid`. Otherwise, we use all cells.
   {if (is.null(targets))
       {message("Reading targets")
        grid = tar_read(pred.grid)
        sat = tar_read(satellite.files)
        model = tar_read(model.full)}
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
