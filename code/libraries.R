suppressPackageStartupMessages(
   {library(stringr)
    library(sf)
    library(ggplot2)
    library(SHAPforxgboost)
    library(Just.universal)
    library(xgboost)
    library(data.table)
    library(pbapply)})
loadNamespace("fst")    # Used for targets with `format = "fst_dt"`
loadNamespace("quarto") # Used by `tarchetypes::tar_quarto`
