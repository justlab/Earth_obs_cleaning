# syntax=docker/dockerfile:1

FROM rocker/r-ver:4.6.0
SHELL ["/bin/bash", "-c"]

RUN apt-get -qq update
RUN apt-get -qq install \
    libgdal-dev libudunits2-dev libssl-dev libglpk-dev libxt-dev \
    curl pandoc libfontconfig1-dev cmake libharfbuzz-dev libfribidi-dev \
    libcairo2-dev libuv1-dev
      # libxt doesn't seem to be strictly necessary, but having it
      # prevents a warning from R's `grSoftVersion()`.
RUN curl >quarto.deb -L https://github.com/quarto-dev/quarto-cli/releases/download/v1.9.37/quarto-1.9.37-linux-amd64.deb && \
    apt-get -qq install ./quarto.deb && \
    rm quarto.deb
RUN apt-get clean

WORKDIR /data/src
COPY code/docker_Rprofile_extra.R ex
RUN cat ex >>$R_HOME/etc/Rprofile.site && rm ex
RUN R -q -e 'install.packages("remotes")'
RUN R -q -e 'remotes::install_github("rstudio/renv@v1.2.2")'

ENV RENV_PATHS_ROOT=/data/renv
ENV OMP_NUM_THREADS=1
  # Work around https://github.com/dmlc/xgboost/issues/2094
ENTRYPOINT ["R", "-q"]
