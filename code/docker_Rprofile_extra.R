options(
    warn = 1,
    repos = c(CRAN = "https://cloud.r-project.org"))

if (interactive())
   {utils::loadhistory(file = "/data/R_history")
    .Last <- function() try(savehistory("/data/R_history"))}
