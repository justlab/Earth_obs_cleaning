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

# These are used by `xgboost.dart.cvtune`.
loadNamespace("ParamHelpers")
loadNamespace("lhs")

# These are used by `download`.
loadNamespace("DBI")
loadNamespace("RSQLite")
loadNamespace("digest")
