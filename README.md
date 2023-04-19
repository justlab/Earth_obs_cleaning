Earth_obs_cleaning is an R-based tool to predict the difference between AERONET and MCD19A2 satellite AOD, and apply that prediction to improve remotely sensed AOD.

# Installation

A Unix-like operating system is assumed. We've developed the project on Ubuntu. A high-end machine may be required for running large workflows, but shouldn't be necessary for the test workflow described in the following section.

- Set some environment variables:
  - `EARTHDATA_USERNAME`, `EARTHDATA_PASSWORD`: [NASA Earthdata login credentials](https://urs.earthdata.nasa.gov/users/new), for downloading the satellite data.
  - `OMP_NUM_THREADS`: This must be set to `1` to work around [an OpenMP bug](https://github.com/dmlc/xgboost/issues/2094).
- Clone the Earth_obs_cleaning repository and `cd` into it.
- Create a configuration file. It should be named `config.yaml` and be placed as a direct child of the repository, alongside this README. See the directory `example` for examples.
  - `data.dir` should be the path to a directory to store downloaded data and cached results. Create it if it doesn't already exist. Expect gigabytes, if not terabytes, of stuff to go in here.
- Start R interactively.
  - `renv::load()`
  - `renv::restore()`
  - `y`
    - This step will install a lot of dependencies. Some system libraries (e.g., `libglpk-dev` on Debian) are required by some of these required R packages. The required packages depend on your OS, and we're not even sure which are required on Ubuntu, so you get to figure that out yourself. Have fun!
  - Quit R, to allow any upgraded packages to be loaded next time.

# The test workflow

Ensure `test.small.daterange` in the configuration file is `TRUE`. Then, in an interactive R session, say `source("code/globals.R")`. Say `tar_make(cv)` to try cross-validation with a few days of data. This is pretty fast, taking only a few minutes, aside from downloading the data. Use `tar_read` (as in `tar_read(cv)`) to see the results.
