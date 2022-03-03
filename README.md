# Earth_obs_cleaning

Predict the difference between AERONET and MCD19A2 satellite AOD, and apply that prediction to improve remotely sensed AOD.

The workflow is contained in `_targets.R`. Rather than using `tar_make()` to run the workflow directly, use a run script such as [R/run_predict.R](R/run_predict.R) which can control which parts of the workflow are best run sequentially or in parallel.
