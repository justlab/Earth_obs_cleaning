suppressPackageStartupMessages(
   {library(stringr)
    library(sf)
    library(ggplot2)
    library(SHAPforxgboost)
    library(Just.universal)
    library(xgboost)
    library(data.table)
    library(pbapply)})
loadNamespace("fst")  # Used for targets with `format = "fst_dt"`
