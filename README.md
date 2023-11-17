Earth_obs_cleaning is an R-based tool to predict the difference between AERONET and MCD19A2 satellite AOD, and apply that prediction to improve remotely sensed AOD.

# Rough notes for the rewrite of this document

1. Put your configuration file in the data directory (the one you're using for the bind mount)

2. `mkdir /tmp` and `/writing` in the data directory.

3. `sudo docker build --tag=earth_obs_cleaning . && sudo docker run --rm -it --mount type=bind,src=/data-coco/Earth_obs_cleaning/docker-data,target=/data -e EARTHDATA_USERNAME -e EARTHDATA_PASSWORD earth_obs_cleaning`

4. `invisible(capture.output({renv::load(); renv::restore(prompt = F)})); source("code/globals.R")`

5. `tar_make(…)`
