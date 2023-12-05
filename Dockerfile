# syntax=docker/dockerfile:1

FROM rocker/r-ver@sha256:fab164fc3015cfeb810a5977c1bef8f8385f75a1f7623eca045b9f4be6f8872a
RUN apt-get -qq update
RUN apt-get -qq install \
    libgdal-dev libudunits2-dev libssl-dev libglpk-dev libxt-dev \
    curl pandoc
      # libxt doesn't seem to be strictly necessary, but having it
      # prevents a warning from R's `grSoftVersion()`.
RUN curl >quarto.deb -L https://github.com/quarto-dev/quarto-cli/releases/download/v1.3.450/quarto-1.3.450-linux-amd64.deb && \
    apt-get -qq install ./quarto.deb && \
    rm quarto.deb
RUN apt-get clean

WORKDIR /src
COPY code/docker_Rprofile_extra.R ex
RUN cat ex >>$R_HOME/etc/Rprofile.site && rm ex
RUN R -q -e 'install.packages("remotes")'
RUN R -q -e 'remotes::install_github("rstudio/renv@v1.0.2")'

RUN R -q -e 'renv::init()' && rm .Rprofile
COPY renv.lock renv.lock
COPY writing writing
COPY code code

ENV RENV_PATHS_ROOT=/data/renv TMPDIR=/data/tmp
ENV OMP_NUM_THREADS=1
  # Work around https://github.com/dmlc/xgboost/issues/2094
ENTRYPOINT ["R", "-q"]
