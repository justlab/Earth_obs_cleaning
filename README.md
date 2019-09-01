# AOD_CONUS
ML AOD CONUS 2019


The drake plan takes __`AOD_CONUS_Report.Rmd`__ and generate html report. I suggest download the repo locally to open the html report. 
```{r}
git clone https://github.com/justlab/AOD_CONUS.git
```
If you are on server coco, [here is the link to the report.](http://coco.5e102.mountsinai.org:8787/files/%2Fhome/liuy29/liuyanguu/AOD_CONUS/AOD_CONUS_Report.html)
  
The __R/__ folder contains all the other code.
  
__`R/drake_plan_conus.Rmd`__ is The drake plan for CONUS;  
__`R/drake_funcs.R`__ contains all the functions used in the drake plans. 
  
The __R/non_drake_code_and_md_files/__ folder contains older code, including `03.2_AOD Variation` about Aeronet stations AOD variation and `plot_predicted_aod470.Rmd` which plots AOD470nm interpolation/exterpolation.


# Notes

Be careful when using data.table with drake. If you modify a data table in a function, this upsteam dataset is modified unless you use data.table::copy to do hard copy in this function. Otherwise this function will rerun again next time you run the workflow. Sometimes drake also gives warnings in such cases.