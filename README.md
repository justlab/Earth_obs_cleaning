Earth_obs_cleaning is an R-based tool to predict the difference between AERONET and MCD19A2 satellite AOD, and apply that prediction to improve remotely sensed AOD.

(The below hasn't been updated since 2b7c939087fb483372f3c339e605515d41ec3388.)

# Installation

A Unix-like operating system is assumed. We've developed the project on Ubuntu. A high-end machine may be required for running large workflows, but shouldn't be necessary for the test workflow described in the following section.

- Set some environment variables:
  - `EARTH_OBS_CLEANING_DATA_DIR`: The path to a directory to store downloaded data and cached results. Create it if it doesn't already exist. Expect gigabytes, if not terabytes, of stuff to go in here.
  - `EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD`: [NASA Earthdata login credentials](https://urs.earthdata.nasa.gov/users/new), for downloading the satellite data.
  - `OMP_NUM_THREADS`: This must be set to `1` to work around [an OpenMP bug](https://github.com/dmlc/xgboost/issues/2094).
  - `EARTH_OBS_CLEANING_NTHREADS` (optional): The number of workers to use for some of the concurrent steps.
- Clone the Earth_obs_cleaning repository and `cd` into it.
- Start `R` interactively.
  - `source("renv/activate.R")`
  - `renv::restore()`
  - `y`
    - This step will install a lot of dependencies. Some system libraries (e.g., `libglpk-dev` on Debian) are required by some of these required R packages. The required packages depend on your OS, and we're not even sure which are required on Ubuntu, so you get to figure that out yourself. Have fun!
  - Quit `R`, to allow any upgraded packages to be loaded next time.

# The test workflow

Set the environment variable `EARTH_OBS_CLEANING_TEST_SMALL_DATERANGE` to `1`. Then, in an interactive R session, say `source("renv/activate.R")`. Say `targets::tar_make(initial_cv_l2_terra_conus)` to try cross-validation with a few days of data, or `targets::tar_make(pred_out_1_test_terra_conus)` to try making new predictions for a single day. These are pretty fast, taking only a few minutes, aside from downloading the data. Use `tar_read` (as in `targets::tar_read(initial_cv_l2_terra_conus)`) to see the results.
