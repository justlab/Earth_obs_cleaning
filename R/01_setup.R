
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check for installed and outdated packages without attaching their namespace ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

package_install <- function(pkg, version, github_source = NULL){
  if(!is.null(github_source)){
    remotes::install_github(github_source)
  } else {
    install.packages(pkg)
}}
check_and_install <- function(pkg, version, github_source = NULL){
  tryCatch(
    if(packageVersion(pkg) < version) package_install(pkg, version, github_source), 
    error = function(e) {
      if(substr(e$message, 1, 26) == "there is no package called"){
        package_install(pkg, version, github_source)
      } else {
        message(e)
      }
})}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries                                  ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(!require(remotes)) install.packages("remotes")
check_and_install("drake", "7.4")
check_and_install("SHAPforxgboost", "0.0.2", "liuyanguu/SHAPforxgboost")

suppressPackageStartupMessages({
  library(sf)
  library(magrittr)
  library(drake)
  library(data.table)
  library(fst)
  library(future)
  library(here)
  library(nngeo)
  library(ggplot2)
  library(SHAPforxgboost)
  library(Just.universal)
})