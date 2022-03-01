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
#' @param absolute whether to use absolute loss (rather than square loss)
#' @param ... other arguments
#'
run.k.fold.cv <- function(k_fold, dataXY_df, y_var,
                          index_train, index_test, xgb_threads, by_var, n_rounds,
                          progress = TRUE, seed = 1234, absolute = F, ...){
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
      objective = (if (absolute) "logcosh" else "reg:squarederror"),
      eval_metric = (if (absolute) "mae" else "rmse"),
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


#' internal function to get threads to use
#'
get.threads <- function(){
  parallel::detectCores()%/%2
}

#' A wrapped function to run xgboost model.
#'
#' @param X predictors, Notice that **should NOT contains Y**.
#' @param Y dependent variable Y.
#' @param xgb_param a list of hyperparameters selected, contains nrounds
#' @return xgboost model object
#'
rfe.fit <- function(X, Y, xgb_param){
  if (!is.null(xgb_param$seed)) set.seed(xgb_param$seed) else set.seed(1234)
  xgb_threads <- get.threads()
  # message(paste("xgb_param is", unlist(xgb_param)))
  xgboost::xgboost(data = as.matrix(X),
                    label = as.matrix(Y),
                    params = xgb_param[names(xgb_param) != 'nrounds'],
                    nrounds = xgb_param$nrounds,
                    verbose = FALSE,
                    nthread = xgb_threads,
                    early_stopping_rounds = 8)
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

# Help func. --------------------------------------------------------------
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
