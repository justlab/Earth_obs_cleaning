# This file has code called by `slurm.sh` and support code for the
# Slurm jobs.

source('code/libraries.R')
source('code/data.R')
source('code/modeling.R')
source('code/util.R')

divisions = \()
  # Group all the cells in the study area into sets ("divisions") of roughly
  # equal size.
   {if (!file.exists(intermediate.path("divisions.parquet")))
       {d = as.data.frame(tar_read(pred.grid)$tile, cells = T, xy = T)
        setDT(d)
        d[, cell := as.integer(cell)]
        d = d[in.sf(x, y, terra::crs(tar_read(pred.grid)), tar_read(region.shape))]
        d[, c("x", "y") := NULL]
        setkey(d, tile, cell)
        d[, division := 1L + ((.I - 1L) %/% 2.5e6L)]
        setkey(d, division, tile, cell)
        setcolorder(d)
        arrow::write_parquet(d, intermediate.path("divisions.parquet"))}
    arrow::read_parquet(intermediate.path("divisions.parquet"))}

make.preds = \(i.work.unit)
   {divisions = arrow::read_parquet(intermediate.path("divisions.parquet"))
    work.unit = `[`(
        CJ(
            year = 2000 : 2025,
            month = 1 : 12,
            division = unique(divisions$division)),
        year > 2000 | month > 2)
            # Terra begins on 2000-02-24
    work.unit = as.list(work.unit[i.work.unit])
    message("Making predictions for work unit ", i.work.unit)
    dput(work.unit)
    dt = with(work.unit, lubridate::make_datetime(year, month))

    d = new.preds.compact(
        dt.start = dt,
        dt.end = lubridate::rollforward(dt, roll_to_first = T),
        cells = divisions[.(work.unit$division), sort(cell)])

    arrow::write_parquet(d, file.path(data.dir, "output",
        with(work.unit, sprintf("%d_%02d_d%d.parquet",
            year, month, division))))}
