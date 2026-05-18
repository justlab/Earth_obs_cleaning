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

Here's how you could interactively make some new predictions for a given range of dates and a given set of points:

```R
source('code/libraries.R')
source('code/data.R')
source('code/modeling.R')
d = new.preds.compact(
    dt.start = as.POSIXct("2024-07-10"),
    dt.end = as.POSIXct("2024-07-16"),
    lonlats = data.frame(
        lon = c(-74.1, -74.1, -74.2),
        lat = c(40.7, 40.8, 40.9)))
```

# VIIRS tiles

VIIRS MAIAC products, referred to as `v19a2` in Earths_obs_cleaning, have a redundant tiling scheme in which multiple tiles can cover the same area. Thus, for `v19a2`, the user must specify the tiles to model in the configuration file, rather than Earths_obs_cleaning choosing tiles automatically to cover the requested region.

To help you choose the tiles you want, here's how you can construct a GeoJSON file for the tile geometry:

```R
library(data.table)
library(sf)

d = fread("tiles.csv")[, -c("V1", "South", "North", "West", "East")]
  # This file can be obtained from https://web.archive.org/web/20260515142741/https://lpdaac.usgs.gov/documents/2436/Zonal_Sinusoidal_Projection_Tile_Definitions.csv
setnames(d, c("h", "v",
    "UL_lat", "UR_lat", "LR_lat", "LL_lat", "UL_lon", "UR_lon", "LR_lon", "LL_lon"))

st_write(dsn = "tiles.geojson", cbind(
    st_as_sf(crs = crs.lonlat, list2DF(list(geometry = apply(simplify = F,
        d[, .(UL_lon, UL_lat, UR_lon, UR_lat, LR_lon, LR_lat, LL_lon, LL_lat, UL_lon, UL_lat)],
        1,
        \(row) st_polygon(list(matrix(row, ncol = 2, byrow = T))))))),
    d[, .(h, v)]))
```

# License

This program is copyright 2019–2025 Kodi B. Arfer, Allan C. Just, Yang Liu, and Johnathan Rush.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [GNU General Public License](http://www.gnu.org/licenses) for more details.
