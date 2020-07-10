
conus_plan_cv <- drake_plan(
  # 5. CV -----------------------------------------------------------
  new_var = target(read_new_var(sat), transform = map(sat = !!sats)),
  new_dt = target(prepare_dt(new_var), transform = map(new_var, .id = FALSE)),
  cv_results = target(run_cv(dt = new_dt), transform = map(new_dt, .id = sat)),
  
  # 6. Report -----------------------------------------------------------
  print1 = target(report_rmse(cv_results),
                  transform = map(cv_results,.id = FALSE)),
  plot1 = target(plot_rmse(cv_results),
                 transform = map(cv_results,.id = FALSE)),
  # figures
  plotAOD = target(plot_AOD(cv_results),
                   transform = map(cv_results,.id = FALSE)),
  plotscatter = target(plot_scatter(cv_results),
                       transform = map(cv_results,.id = FALSE)),
  
  shap_long = target(get_SHAP_long(cv_results),
                     transform = map(cv_results,.id = FALSE)),
  plotSHAP1 = target(plot_SHAP(shap_long),
                     transform = map(shap_long,.id = FALSE)),
  plotSHAP2 = target(plot_SHAP_scatter(cv_results, shap_long),
                     transform = map(cv_results, shap_long,.id = FALSE)),
  report = rmarkdown::render(
    knitr_in(here("AOD_CONUS_Report.Rmd")),
    output_file = file_out(here("AOD_CONUS_Report.html")),
    output_format = rmarkdown::html_document(toc = TRUE))
)
conus_cv_config <- drake_config(conus_plan_cv, cache = ssd_cache)
