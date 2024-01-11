options(
    warn = 1,
    repos = c(CRAN = "https://cloud.r-project.org"))

if (interactive())
   {utils::loadhistory(file = "/data/R_history")
    .Last <- function() try(savehistory("/data/R_history"))}

if (file.exists("/data/R-packages-installed"))
   {invisible(utils::capture.output({renv::load(); renv::restore(prompt = F)}))
    source("code/globals.R")}
