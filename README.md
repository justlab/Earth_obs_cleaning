Earth_obs_cleaning is an R- and Docker-based tool to predict the difference between AERONET and MCD19A2 satellite aerosol optical depth (AOD), and apply such predictions to improve the same AOD product.

# Installation

A high-end machine may be required for running large workflows, but shouldn't be necessary for the test workflow described in the following section.

1. Create a data directory. The data directory will store code, downloaded data, cached results, and temporary files. For large workflows, expect terabytes of stuff to go into it.
2. Create a configuration file named `config.yaml` in the data directory. See the directory `example` in the Earth_obs_cleaning repository for examples.
3. Clone this repository to the data directory, and name the new directory `src`. (Actually, of the items in this repository, only `renv.lock`, `code`, and `writing` are required.)
4. To build the Docker image, use the command `docker build --tag=earth_obs_cleaning .`
5. To use the image to create a container and start R interactively, say `docker run --rm -it --mount type=bind,src=DPATH,target=/data -e EARTHDATA_USERNAME -e EARTHDATA_PASSWORD earth_obs_cleaning`
    - Replace `DPATH` with the path to your data directory.
    - Notice that the environment variables `EARTHDATA_USERNAME` and `EARTHDATA_PASSWORD` should be set in your real environment; these are [NASA Earthdata login credentials](https://urs.earthdata.nasa.gov/users/new) for downloading satellite data.
    - `--rm` is used to automatically delete the container after the R process exits. This is convenient but not necessary.
6. In R, say:
   - `renv::init()`
     - `1`
   - `unlink(".Rprofile")`
   - `cat("TRUE\n", file = "/data/R-packages-installed")`

# Usage

Run a Docker container as described above in step 5 above. (If you installed `renv` packages in this R session, quit and restart.) You can now use `tar_make` to build targets.

To run the test workflow, ensure `test.small.daterange` in the configuration file is `TRUE`. Then say `tar_make(cv)` to try cross-validation with a few days of data. This is pretty fast, taking only a few minutes, aside from downloading the data. Use `tar_read` (as in `tar_read(cv)`) to see the results.
