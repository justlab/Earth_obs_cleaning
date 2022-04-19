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
