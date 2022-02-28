# Contains most functions required for the AOD cleaning targets workflow

#' Return AERONET points that intersect the polygon describing the area of
#' interest, and include the cell ID (\code{idM21pair0}) from the reference
#' CONUS MODIS grid.
#'
#' @param stations data.table including coordiantes (\code{lon, lat}) points of
#'   AERONET stations
#' @param reg_polygon SF polygon defining the area of interest
#' @param refgrid_path FST with coordinates of all raster cells in the area of
#'   interest
#' @param refras_path path to MODIS reference raster stack
#' @return a subset of geometry points within the given polygon
select_stations <- function(stations, reg_polygon, refgrid_path, refras_path){
  aerpts = st_as_sf(stations, coords = c("lon", "lat"), crs = 4326)
  if(st_crs(aerpts) != st_crs(reg_polygon)){ # polygons are in nad83
    aerpts <- st_transform(aerpts, crs = st_crs(reg_polygon))
  }
  aerpts <- aerpts[st_intersects(aerpts, reg_polygon, sparse = FALSE),]
  aerpts <- station_cell_ids(aerpts, refgrid_path, refras_path) # now in sinusoidal CRS
  aerpts
}

#' Get the unique ID of the grid cell an AERONET station is in.
#'
#' @param stations_sf SF points of AERONET stations
#' @param refgrid_path FST with coordinates of all raster cells in the area of
#'   interest
#' @param refras_path path to MODIS reference raster stack
#' @return SF points of AERONET stations with unique MODIS grid cell IDs added
#'   (\code{idM21pair0})
station_cell_ids = function(stations_sf, refgrid_path, refras_path){
  gras = raster(refras_path, band = 4)
  stations_sf = st_transform(stations_sf, crs = st_crs(gras))
  stations_sf$cell_index = raster::extract(gras, as(stations_sf, 'Spatial'))

  gDT = read_fst(refgrid_path, as.data.table = TRUE, columns = c('idM21pair0', 'cell_index'))
  setkey(gDT, cell_index)
  setDT(stations_sf)
  setkey(stations_sf, cell_index)
  stations_sf[gDT, idM21pair0 := idM21pair0]
  stations_sf[, cell_index := NULL]
  setkey(stations_sf, Site_Name)
  stations_sf = st_sf(stations_sf)
  # previous row names no longer meaningful after resorting
  attributes(stations_sf)$row.names <- 1:nrow(stations_sf)
  stations_sf
}

#' Load Aeronet AOD measurement data from text files sourced from downloaded
#' \code{tar.gz} file. Only load the data for the specified stations.
#'
#' @param aod_dir directory where the \code{.lev20} files were extracted to.
#' @param stations SF points of AERONET stations
#' @param date_start subset to observations this date or later
#' @param date_end subset to observations this date or earlier
#' @return data.table of AERONET observations, joined with MODIS reference grid
#'   unique ID
get_stn_data <- function(aod_dir, stations, date_start = NULL, date_end = NULL){
  stn_names = unique(stations$Site_Name)
  # open files containing names of stations in the specified region
  aer_files_dir <- sapply(paste0(unique(stn_names),".*\\.lev20"), FUN = list.files,
                          path = aod_dir, full.names = TRUE)
  found_files <- sapply(aer_files_dir, function(x) length(x) > 0)
  aer_files_dir = aer_files_dir[found_files]
  if(length(aer_files_dir) > 0){
    t0 <- fread(aer_files_dir[[1]], nrows = 10) # to get variable names:
    vars_aod <- intersect(grep("AOD_", names(t0), value = T), grep("nm", names(t0), value = T))
    # sort the wave lengths varnames from low to high, and update the vars_aod
    x_nm <- sort(readr::parse_number(vars_aod)) # 340, 380, ... , 1640
    vars_aod <- paste0("AOD_", x_nm,"nm") # the AOD_...nm variables
    vars_wv <- paste0("Exact_Wavelengths_of_AOD(um)_", x_nm,"nm") # the Exact wv length

    # data_path is a single file path
    read_aod <- function(data_path, date_start = NULL, date_end = NULL){
      dt <- fread(data_path, select =  c(vars0, vars_aod, vars_wv))
      dt[, stn_time := as.POSIXct(paste(`Date(dd:mm:yyyy)`, `Time(hh:mm:ss)`),
                                  format = "%d:%m:%Y %H:%M:%S", tz = "UTC")]
      dt[, c("Date(dd:mm:yyyy)", "Time(hh:mm:ss)") := NULL]
      for (i in names(dt)) dt[get(i) == -999, (i):= NA] # set NA
      if(!is.null(date_start)) {dt = dt[as.Date(stn_time) >= as.Date(date_start), ] }
      if(!is.null(date_end))   {dt = dt[as.Date(stn_time) <= as.Date(date_end), ] }
      setcolorder(dt, c("AERONET_Site_Name", "stn_time"))
      dt
    }

    file_list <- lapply(unlist(aer_files_dir), read_aod, date_start = date_start, date_end = date_end)
    aer_data <- rbindlist(file_list)

    aer_data[, aer_date := as.Date(stn_time)]
    aer_data
  } else {
    NULL
  }
}

#' Generate a tibble with a column for year and a column with a list of every
#' date in the year
#'
#' @param years vector of four-digit years
#' @return a tibble with a "year" column and a "dates" column containing a list
#'   of every date in the year
dates_year <- function(years){
  dates = lapply(years, function(y) seq.Date(as.Date(paste0(y, '-01-01')),
                                            as.Date(paste0(y, '-12-31')), 1))
  tibble::tibble(year = years, dates = dates)
}

dates_year_list <- function(years){
  dates = lapply(years, function(y) seq.Date(as.Date(paste0(y, '-01-01')),
                                            as.Date(paste0(y, '-12-31')), 1))
  names(dates) <- years
  dates
}

#' Assign a month index from 1970-01-01 to all observation dates.
#'
#' Month indexes are used to batch the extraction of training data into chunks and increase targets throughput.
#' @param aer_data AERONET observation data
#' @return data.table of AERONET observation data with a monthid column added
assign_monthid <- function(aer_data){
  start = as.Date('1970-01-01') - 1
  aer_data[, monthid := as.period(as.Date(aer_date) - start) %/% months(1)]
  aer_data
}

#' Filter data.table of AERONET observations to only those in the given dates
#'
#' @param aer_data AERONET observation data.table with a date column `aer_date`
#' @param dates vector of dates to filter the AERONET observation
#' @return data.table of AERONET observations in the given dates
#'
filter_aer_bydate <- function(aer_data, dates){
  aer_data[aer_date %in% dates, ]
}

#' prepare AERONET observation data for interpolation of AOD470nm
#'
#' interpolating to wavelength 470 nm, input as 0.47
#'
#' @param aer_data the dataset contains all the AOD measurements and exact wavelengths
#' @param aer_stns simple feature collection of station locations
#' @return dataset of aer_data with pred, also with `nearest_refgrid` next to Site_Name
interpolate_aod <- function(aer_data, aer_stns){
  keepcols = grep('1020|1640', names(aer_data), invert = TRUE)
  aer_data <- aer_data[, ..keepcols]
  aer_data[, obs:= .I]
  setkey(aer_data, obs)
  vars_wv <- grep("Exact_Wavelengths_of_AOD", names(aer_data), value = TRUE)
  vars_aod_sub <- grep("AOD_", names(aer_data), value = TRUE)
  aer_data$N_NAAOD <- rowSums(!is.na(aer_data[, ..vars_aod_sub]))
  # create long-format data
  d <- melt(aer_data[, c(vars_aod_sub, vars_wv, "obs"), with = FALSE],
            measure = list(vars_aod_sub, vars_wv), value.name = c("aod", "exact_wv"))
  # choose non-NA and >0 AOD
  d <- d[!is.na(aod) & aod > 0]
  # ids for wavelengths observed on one side of 470nm
  obs_os <- d[, by = obs, .(oneside = (max(exact_wv)-0.47)*(min(exact_wv)-0.47)>0)][oneside == TRUE, obs]
  # the predict part:
  d2 <- d[, by = obs,
     {m = lm(log(aod) ~ log(exact_wv) + I((log(exact_wv))^2))
      if (m$rank == length(coef(m)))
        # Only use the model if it's full-rank.
          list(AOD_470nm = exp(predict(m,
              data.frame(exact_wv = 0.47))))}] # notice to put 0.47 not 470
  setkey(d2, obs)
  aer_data_wPred <- d2[aer_data]
  # set NA: At least 4 observations wanted for exterpolation
  aer_data_wPred[N_NAAOD < 4 & obs %in% obs_os, AOD_470nm := NA]
  # remove unnecessary variables:
  aer_data_wPred <- aer_data_wPred[, -c(vars_aod_sub, vars_wv), with = FALSE]

  # join station locations to observations
  setnames(aer_data_wPred, "AERONET_Site_Name", "Site_Name")
  setkey(aer_data_wPred, Site_Name)
  aer_sites = as.data.table(aer_stns)
  aer_sites[, c('x_sinu', 'y_sinu') := as.data.table(st_coordinates(geometry))]
  aer_sites = aer_sites[, .(Site_Name, x_sinu, y_sinu)]
  setkey(aer_sites, Site_Name)
  aer_data_wPred <- aer_sites[aer_data_wPred]
  setkey(aer_data_wPred, Site_Name, stn_time)
  aer_data_wPred
}

#' Return the SF points for stations reporting data in the previously-selected
#' time period
#'
#' @param stations SF points of AERONET stations
#' @param stn_data data.table of station observations for previously-selected time period
filter_stations = function(stations, stn_data){
  stations[stations$Site_Name %in% stn_data$Site_Name, ]
}

#' Get MODIS grid unique IDs (\code{idM21pair0}) for every cell within specified
#' distance from AERONET stations.
#'
#' @param stations
#' @param refgrid_path FST with coordinates of all raster cells in the area of interest
#' @param dist_km Return cell IDs within this distance, in kilometers
cells_in_buffer = function(stations, refgrid_path, dist_km = 270){
  # set up so that function can either receive a single station (1 row) or
  # multiple stations, allowing external parallelism from dynamic branching
  gDT = read_fst(refgrid_path, as.data.table = TRUE, columns = c('idM21pair0', 'x_sinu', 'y_sinu'))

  cells_by_station = function(rownum, gDT){
    stn_coords = st_coordinates(stations[rownum, ])
    subDT = gDT[x_sinu <= stn_coords[, 'X'] + dist_km * 1000 &
                x_sinu >= stn_coords[, 'X'] - dist_km * 1000 &
                y_sinu <= stn_coords[, 'Y'] + dist_km * 1000 &
                y_sinu >= stn_coords[, 'X'] - dist_km * 1000]
    subDT[, Site_Name := stations[rownum, ]$Site_Name]
    subDT[, site_dist := as.numeric(st_distance(st_as_sf(.SD[, .(x_sinu, y_sinu)],
                                                         coords = c(1,2), crs = st_crs(stations)),
                                     stations[rownum, ])[, 1])]
    subDT[site_dist <= dist_km * 1000, .(Site_Name, idM21pair0, site_dist)]
  }
  # would exceed row maximum after ~8000 stations with 270km distance
  rbindlist(lapply(1:nrow(stations), FUN = cells_by_station, gDT = gDT))
}

#' Return a reference raster with the extent of the FST reference grid instead of
#' the reference raster TIF.
#' This will be compatiable with `aoiname = 'conus'`, but would require further
#' cropping to use with other AOIs.
#'
crop_refras_mcd <- function(refgrid_path, mcd19path,
                            ref_uid = 'idM21pair0', aoiname = 'conus'){
  if(!aoiname %in% c('conus', 'nemia')) stop('Only CONUS region has been implemented')
  # get the mcd19 file for the first date in an arbitrary year
  mcdDT = read_fst(list.files(file.path(mcd19path, 2010), '*.fst',
                              full.names = TRUE)[1],
                   as.data.table = TRUE, columns = ref_uid)
  rg = read_fst(refgrid_path, as.data.table = TRUE,
                columns = c(ref_uid, 'x_sinu', 'y_sinu', 'cell_index'))
  setkeyv(mcdDT, ref_uid)
  rasterFromXYZ(rg[mcdDT, c('x_sinu', 'y_sinu', 'cell_index'), with = FALSE],
                crs = crs_sinu)
}

#' Return a matrix classifying cells as being within a distance from the center.
#'
#' Assumes square matrix with odd number of rows & columns. The width (and
#' height) will be 2x the radius plus one cell.
#'
#' @param radius distance, in kilometers, from center cell to classify as inside
#'   circle
#' @param cellsize dimensions of raster cells in meters. Pixels assumed square.
#' @param matrix if TRUE, return a logical matrix. If FALSE (default), return a
#'   data.table with two columns.
#' @return a data table with offsets from the central cell that are within the
#'   specified radius
circle_mat = function(radius, cellsize = 926.6254, matrix = FALSE){
  radius = radius * 1000
  width = ceiling(radius/cellsize*2+1)
  if(width%%2 == 0) width <- width + 1 # ensure odd row and column count
  height = width # assuming square matrix
  circle_df = expand.grid(mrow = 1:height, mcol = 1:width)
  center_val = median(circle_df$mrow) # assuming square matrix
  circle_df$dist_cells = sqrt((circle_df$mrow - center_val)^2 + (circle_df$mcol - center_val)^2)
  circle_df$circ = ifelse(circle_df$dist_cells * cellsize <= radius, TRUE, FALSE)

  if(matrix == TRUE){
    circle_mat = matrix(as.numeric(circle_df$circ), height, width)
    circle_mat
  } else {
    setDT(circle_df)
    circle_df[, offset_x := mcol - center_val]
    circle_df[, offset_y := mrow - center_val]
    circle_df[circ == TRUE, .(offset_x, offset_y)]
    circle_df
  }
}

# Extent to create a raster version of focal filter matrix
# May extend beyond the AOD raster
get_focal_extent = function(x, c1, r1, radw){
  center_x = xFromCol(x, c1)
  center_y = yFromRow(x, r1)
  r = res(x)
  xn <- center_x - ((0.5 * r[1]) + (r[1] * radw))
  xx <- center_x + ((0.5 * r[1]) + (r[1] * radw))
  yn <- center_y - ((0.5 * r[2]) + (r[2] * radw))
  yx <- center_y + ((0.5 * r[2]) + (r[2] * radw))
  ext(c(sort(c(xn, xx))), sort(c(yn, yx)))
}

#' Roll join a single day of AERONET data to nearest MCD19A2 overpass and
#' calculate derived values from MCD19A2 within specified distances from the
#' AERONET station.
#'
#' @param aer_data A single day of AERONET observations
#' @param load_sat input "terra" or "aqua"
#' @param buffers_km vector of buffer radii in kilometers
#' @param aer_stn data table of AERONET station names and reference unique ID
#'   for satellite AOD cells
#' @param hdf_root path to MCD19A2 HDF files
#' @param agg_level how much to aggregate the input MODIS AOD to use for large
#'   focal radii
#' @param agg_thresh use the aggregated AOD when `(radius * 1000 / agg_level / mcd_res) > agg_thresh`,
#'   where `radius` is in km, and `mcd_res` (resolution of input AOD) is
#'   in meters.
#' @param vrt_path directory to store overpass VRTs that reference HDF files.
#' @param rolldiff_limit maximum time between the AERONET and satellite
#'   observations
derive_mcd19_vars = function(aer_data, load_sat, buffers_km, aer_stn, hdf_root,
                             agg_level, agg_thresh, vrt_path,
                             rolldiff_limit = as.difftime(30, units = 'mins')){
  # 1. Prepare AERONET data
  aer_stn = st_sf(aer_stn) # testing whether passing in a DT version avoids error with vctrs package
  sv_aer = vect(aer_stn)
  setDT(aer_data)
  if(nrow(aer_data) == 0){
    return(data.table(NA))
  }
  aer_data = interpolate_aod(aer_data, aer_stn)

  # Get the date of this chunk of AERONET data, find HDFs from same date
  if(aer_data[, uniqueN(aer_date)] > 1) stop('Use tar_group_by and pattern=map to send one date at a time to derive_mcd19_vars')
  this_date = aer_data[1, aer_date]

  hdf_paths = list.files(file.path(hdf_root, format(this_date, '%Y.%m.%d')),
                         pattern = '\\.hdf$', full.names = TRUE)

  if(length(hdf_paths) == 0){
    return(data.table(NA))
  }
  # 2. Join AERONET to satellite AOD
  binDT = bin_overpasses(hdf_paths)
  binDT = binDT[sat == load_sat]
  if(nrow(binDT) == 0){
    return(data.table(NA))
  }

  day_op = get_overpasses_vrts(hdf_paths, binDT, load_sat, vrt_path)

  # lapply by overpass, which are the groups in day_op
  tryCatch({
  day_rasters = mapply(FUN = rasterize_vrts,
                       op_vrts = day_op,
                       opDT = lapply(1:length(day_op), function(i) binDT[overpass_bin == i]),
                       SIMPLIFY = FALSE)
  }, error = function(e) {
    emsg = paste0('this_date = ', this_date, '\n',
           'length(day_op) = ', length(day_op), '\n',
           'binDT[, uniqueN(overpass_bin)] = ', binDT[, uniqueN(overpass_bin)], '\n',
           'error = ', e)
    stop(emsg)
  })
  rm(day_op)

  # roll join will update value of the time column in X to the value of the time
  # column in i. Copy AERONET's stn_time to a column with the same name as
  # MCD19's time column so we can compare the time difference afterwards using
  # meaningful column names.
  setnames(aer_data, 'stn_time', 'overpass_time')
  aer_data[, stn_time := overpass_time]

  # extract AOD overpass values by each overpass
  overpass_stats <- function(single_op){

    # extract values from each raster for the overpass
    op_vals = lapply(single_op, function(op){
      # only record cell ID for the raster stack with AOD
      if(length(grep('Optical_Depth', names(op))) > 0){
        dt = setDT(terra::extract(op, sv_aer, cells = TRUE))
      } else {
        dt = setDT(terra::extract(op, sv_aer, cells = FALSE))
      }
      setkey(dt, ID)
    })
    op_vals = Reduce(merge.data.table, op_vals)
    op_vals[, Site_Name := aer_stn$Site_Name]
    setnames(op_vals, 'Optical_Depth_047', 'MCD19_AOD_470nm')
    op_vals = op_vals[!is.na(MCD19_AOD_470nm)]
    op_vals[, overpass_time := as.POSIXct(overpass_time, tz = 'UTC', origin = '1970-01-01')]

    # roll join AERONET to satellite AOD
    setkey(op_vals, Site_Name, overpass_time)
    rj <- aer_data[op_vals, roll = 'nearest', nomatch = 0] # using keys to join

    # Only keep roll joins within specified time difference limit
    rj[, rj_difftime := overpass_time - stn_time]
    rj = rj[abs(rj_difftime) <= rolldiff_limit, ]

    if(nrow(rj) > 0){
      # 3. Derived values from focal statistics on satellite AOD

      # find the AOD layer in the list of rasters for the overpass
      mcd_aod = single_op[unlist(lapply(single_op,
        function(s) 'Optical_Depth_047' %in% names(s)))][[1]][['Optical_Depth_047']]
      mcd_nna = !is.na(mcd_aod)
      mcd_res = res(mcd_aod)[1] # assuming square pixels
      # Aggregate AOD raster for use in larger focal radii
      mcd_agg_aod = terra::aggregate(mcd_aod, fact = agg_level, fun = 'mean', na.rm = TRUE)
      mcd_agg_nna = terra::aggregate(mcd_nna, fact = agg_level, fun = 'mean', na.rm = TRUE)

      filter_matrices = list()
      use_agg = list()
      for(radius in buffers_km){
        if(radius * 1000 / agg_level / mcd_res < agg_thresh){
          # use original resolution
          use_agg[[paste0('r', radius)]] <- FALSE
          this_mat = circle_mat(radius, cellsize = mcd_res, matrix = TRUE)
        } else {
          # use aggregated resolution
          use_agg[[paste0('r', radius)]] <- TRUE
          this_mat = circle_mat(radius, cellsize = mcd_res * agg_level, matrix = TRUE)
        }
        this_mat[this_mat == 0] <- NA
        filter_matrices[[paste0('r', radius)]] <- this_mat
      }
      rm(radius)

      # Calculate focal statistics for each unique cell containing an AERONET station
      # returns a data.table to merge (send unique values in case any stations in same cell)
      station_focal <- function(cellid){
        agg_cellid = cellFromXY(mcd_agg_aod, xyFromCell(mcd_aod, cellid))
        focal_list = list(cellid = cellid, agg_cellid = agg_cellid)
        # buffers should be in km, resolution should be in m
        for(radius in buffers_km){
          filter_mat = filter_matrices[[paste0('r', radius)]]
          focw = ncol(filter_mat)  # width of focal matrix (assumed square matrix)
          radw = (focw - 1) / 2    # width of the radius in cells, without central cell
          # Extract matrix of AOD values around station
          if(use_agg[[paste0('r', radius)]] == FALSE){
            # original AOD raster resolution
            crid = rowFromCell(mcd_aod, cellid) # central row id
            ccid = colFromCell(mcd_aod, cellid) # central column id
            filter_raster = rast(filter_mat, crs = crs_sinu,
                                 extent = get_focal_extent(mcd_aod, ccid, crid, radw))
            aod_extract = crop(mcd_aod, filter_raster)
            nna_extract = crop(mcd_nna, filter_raster)
          } else {
            # aggregated AOD raster
            crid = rowFromCell(mcd_agg_aod, agg_cellid) # central row id
            ccid = colFromCell(mcd_agg_aod, agg_cellid) # central column id
            filter_raster = rast(filter_mat, crs = crs_sinu,
                                 extent = get_focal_extent(mcd_agg_aod, ccid, crid, radw))
            aod_extract = crop(mcd_agg_aod, filter_raster)
            nna_extract = crop(mcd_agg_nna, filter_raster)
          }
          filter_crop = crop(filter_raster, aod_extract)
          mult_aod = aod_extract * filter_crop
          mult_nna = nna_extract * filter_crop
          focal_list[[paste0('Mean_AOD', radius, 'km')]] <- mean(mult_aod[], na.rm = TRUE)
          focal_list[[paste0('pNonNAAOD', radius, 'km')]] <- mean(mult_nna[], na.rm = TRUE)
        }
        setDT(focal_list)
        focal_list
      }
      focal_stats = rbindlist(lapply(rj[, unique(cell)], station_focal))
      setkey(focal_stats, cellid)
      setkey(rj, cell)
      rj = rj[focal_stats]

      # difference between central cell satellite AOD and mean AOD in buffers
      for(distx in buffers_km){
        rj[, c(paste0('diff_AOD', distx, 'km')) := MCD19_AOD_470nm - get(paste0('Mean_AOD', distx, 'km'))]
      }
      rj[, c('ID', 'x_sinu', 'y_sinu') := NULL] # remove AERONET station row index and coordinates
      rj
    } else {
      # This will create an extra "V1" column in rowbound final target with all NAs,
      # but it gets around the error of writing NULL to FST format
      data.table(NA)
    }
  }
  outDT = rbindlist(lapply(day_rasters, overpass_stats), fill = TRUE)
  if(ncol(outDT) > 1){
    outDT = outDT[!is.na(Site_Name)]
    if('V1' %in% names(outDT)) outDT[, V1 := NULL]
    outDT
  } else { # no RJ results for any overpass on this day
    data.table(NA)
  }
}

#' Open one day of MCD19A2 data from a FST file. Adds a column for the
#' satellite.
#'
#' @param sat input "terra" or "aqua", if sat !="terra", load aqua
#' @param daynum day number as character with padding to three digits, leading
#'   zeroes
#' @param filepath path to the year of MCD19A2 daily best overpasses as FST
#'   files
#' @return one day of MCD19A2 AOD data as a data.table, keyed by
#'   \code{idM21pair0, overpass_time}
read_mcd19_one <- function(sat = "terra", daynum, filepath, load_year,
                           ref_uid = 'idM21pair0',
                           columns = c("overpass_index", "Optical_Depth_047",
                                       "AOD_Uncertainty", "Column_WV", "AOD_QA",
                                       "RelAZ", "idM21pair0", "overpass_time")){
  columns = unique(c(ref_uid, 'overpass_time', columns))
  if (sat != "terra") choose = "A" else choose = "T" # load terra by default
  mcd_file = file.path(filepath, paste0('mcd19_conus_', choose, '_', load_year,
                                        daynum, '.fst'))
  if(file.exists(mcd_file)){
    dt = read.fst(mcd_file, as.data.table = TRUE, columns = columns)
    dt[, sat := choose]
    #setnames(dt, "idM21pair0", "nearest_refgrid") # rename for join later
    #setkey(dt, nearest_refgrid, join_time)
    if('overpass_time' %in% columns){
      setkeyv(dt, c(ref_uid, 'overpass_time'))
    } else {
      setkeyv(dt, ref_uid)
    }
    dt
  } else {
    NULL
  }
}

create_qc_vars <- function(dt){
  dt[, qa_bits := bitwAnd(AOD_QA, strtoi("111100000000", base = 2))]
  dt[, qa_best := 0]
  dt[qa_bits == 0 , qa_best := 1]

  dt[, qa_bits := bitwAnd(AOD_QA, strtoi("11000", base = 2))]
  dt[qa_bits==0, qa_lwsi := "land"]
  dt[qa_bits==8, qa_lwsi := "water"]
  dt[qa_bits==16, qa_lwsi := "snow"]
  dt[qa_bits==24, qa_lwsi := "ice"]
  dt[, qa_bits:= NULL]
  return(dt)
}

#' Calculate DV, remove rows without DV, and calculate dayint and QC var IVs.
#'
#' Optionally subset the MCD19A2 observations to those with dates in the
#' \code{date_range} vector.
#'
#' @param dt data.table of MCD19A2 observations
#' @param date_range vector of Date type containing all dates to retain in the
#'   subset
#' @return data.table with DV and additional IV columns ready for running CV
prepare_dt <- function(dt, date_range = NULL){
  setnames(dt, "Optical_Depth_047", "MCD19_AOD_470nm", skip_absent=TRUE)
  if(!is.null(date_range)){
    if(class(date_range) == 'list') date_range = date_range[[1]]
    dt = dt[aer_date %in% date_range, ]
  }
  # The dependent variable: diff_AOD = MCD19 - AERONET = Optical_Depth_047 - AOD_470nm
  dt[, diff_AOD := MCD19_AOD_470nm - AOD_470nm]
  dt <- create_qc_vars(dt)
  dt[, dayint := as.integer(as.Date(overpass_time))]
  dt = dt[!is.na(get(y_var))]
  dt
}

#' Do initial CV with hyperparameter selection with DART and rank features.
#'
#' Must provide either stn_var or day_var for binning folds.
#'
#' @param data data.frame of observations with DV and IVs
#' @param y_var column name to predict
#' @param features character vector of column names to train on
#' @param stn_var if binning by stations, provide the station column name
#' @param day_var if binning by days, provide the day column name
# Depends on functions in xgboost_cv_RFE.R
initial_cv_dart <- function(
  data,
  y_var,
  features,
  k_fold = 5,
  n_rounds = 100,
  stn_var = NULL,
  day_var = NULL,
  progress = TRUE
){
  xgb_threads <- get.threads()

  mDT <- data.table::copy(setDT(data)) # needed for drake, uncertain about targets
  if(!is.null(day_var)) {
    mDT[, dayint := as.integer(get(day_var))]
    by_var = "day"
  }
  if(!is.null(stn_var)) {
    mDT[, stn := get(stn_var)]
    by_var = "stn"
  }
  if (!is.null(day_var) & !is.null(stn_var))
    stop("Please provide either day_var or stn_var.")

  # some internal variables names
  y_formula <- as.formula(paste(y_var, "~."))
  y_var_pred <- paste0(y_var, "_pred") # name of the predicted y
  y_var_pred_whole <- paste0(y_var, "_pred_whole") # name of the predicted y

  # bin data into specified number of folds
  temp <- prepare.bin(mDT, by_var = by_var, k_fold = k_fold)
  mDT <- temp$data
  bin_list <- temp$bin # list of values of selected variable (stn or day) by fold
  rm(temp)
  index.fs.list <-  list(
    # return the observations in each fold
    stn = function(x, bin_list) mDT[stn%in%bin_list[[x]], which = TRUE],
    day = function(x, bin_list) mDT[dayint%in%bin_list[[x]], which = TRUE]
  )

  index_train <- index_test <- list()
  for (k in 1:k_fold){
    index_test[[k]]  <- index.fs.list[[by_var]](k, bin_list)
    index_train[[k]] <- -index_test[[k]]
  }

  # run k-fold cv and record SHAP matrix, predicted value, overall rmse...
  if (!y_var%in%features) features <- c(features, y_var)

  message("Run k-fold cv \n")
  cv_results <- run.k.fold.cv(k_fold = k_fold,
                              run_param_cv = TRUE, # run DART to select hyperparameters
                              dataXY_df = mDT[, ..features],
                              by_var = by_var,
                              n_rounds = n_rounds,
                              y_var = y_var,
                              progress = progress,
                              index_train = index_train,
                              index_test = index_test,
                              xgb_threads = xgb_threads)

  mDT <- cbind(mDT, cv_results$y_pred_dt)

  # output
  c(list(by_var = by_var, bin_list = bin_list, mDT_wPred = mDT), cv_results)
}

#' Summarize CV statistics on all initial_cv objects
#'
#' @param cv_list list of all CV output objects; may include every year and both
#'   satellites.
#' @return data.table summarizing CV statistics for each year and satellite
cv_summary <- function(cv_list){
  cv_list = unlist(cv_list, recursive = FALSE)
  stats_list = vector(mode = "list", length = length(cv_list))
  difftimes_list = vector(mode = "list", length = length(cv_list))
  for(i in 1:length(cv_list)){
    cv = cv_list[[i]]
    stats = cv_reporting(cv)
    stats$sat <- str_extract(names(cv_list)[[i]], 'terra|aqua')
    stats$year <- cv$mDT_wPred[1, year(aer_date)]
    stats_list[[i]] <- stats
    difftimes_list[[i]] <- round(summary(as.numeric(abs(cv$mDT_wPred$rj_difftime))),0)
  }
  # prediction summary stats
  statsDT = rbindlist(lapply(stats_list, as.data.table))
  statsDT[, MAE_pct_change :=
          paste0(round((MAE_uncorr-MAE_corr)/MAE_uncorr, 2) * 100, '%')]
  statsDT[, MAD_pct_change :=
          paste0(round((MAD_mcd19-MAD_aodhat)/MAD_mcd19, 2) * 100, '%')]
  setcolorder(statsDT, c('sat', 'year',
                       'MAE_uncorr', 'MAE_corr', 'MAE_pct_change',
                       'rmse', 'MAD_mcd19', 'MAD_aodhat', 'MAD_pct_change'))
  # difftime distribution
  difftimeDT = rbindlist(lapply(difftimes_list, function(x) as.list(x)))
  difftimeDT[, c('sat', 'year') := statsDT[, .(sat, year)]]
  setcolorder(difftimeDT, c('sat', 'year'))

  setkey(statsDT, sat, year)
  setkey(difftimeDT, sat, year)
  list(stats = statsDT, difftimes = difftimeDT)
}

#' Calculate CV statistics on a single \code{initial_cv_dart()} output list object
#'
cv_reporting <- function(cv){
  dt = cv$mDT_wPred
  mae = function(v1, v2) mean(abs(v1 - v2))
  mad = function(v1) mean(abs(v1 - median(v1)))
  dt[, aod_hat := MCD19_AOD_470nm - diff_AOD_pred]

  list(
    MAE_uncorr = round(dt[, mae(MCD19_AOD_470nm, AOD_470nm)],3),
    MAE_corr   = round(dt[, mae(aod_hat, AOD_470nm)],3),
    rmse = round(cv$rmse_all_folds,3),
    MAD_mcd19  = round(mad(dt$MCD19_AOD_470nm),3),
    MAD_aodhat = round(mad(dt$aod_hat),3),
    stn_count = dt[, uniqueN(stn)],
    train_N = dt[, .N]
  )
}

# Prediction ####

#' Create a data.table with every prediction date and the corresponding trained
#' XGBoost model to use.
#'
model_files_by_date <- function(mft, pred_dates){
  mft = mft[, .(year, file_path)]
  dt = data.table(pred_date = pred_dates)
  dt[, year := year(pred_date)]
  outdt = mft[dt, on = 'year', nomatch = 0]

  missing_dates = !pred_dates %in% outdt$pred_date
  if(any(missing_dates)){
    missing_years = paste(unique(year(pred_dates[missing_dates])), collapse = ', ')
    stop('Requested predictions in ', missing_years,
         ' but no model has been trained for those years.')
  } else {
    outdt
  }
}

#' #' Given an input SF object, return a terra:::SpatExtent in a new CRS,
#' #' optionally extended in all directions by `expand_distance`, using the units
#' #' of the `new_crs`.
#' get_aoi_ext <- function(shape, new_crs, expand_distance = 0){
#'   shape = st_transform(shape, new_crs)
#'   bbox = st_bbox(shape)
#'   if(expand_distance != 0){
#'     mod_box = st_bbox(c(xmin = -expand_distance, xmax = expand_distance,
#'                         ymin = -expand_distance, ymax = expand_distance))
#'     bbox = bbox + mod_box
#'   }
#'   ext(bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax)
#' }

#' Given an input SF object, return a terra:::SpatExtent. Optionally reproject
#' to `to_crs`.
sf_to_ext <- function(sf, to_crs = NULL){
  if(!is.null(to_crs)) sf = st_transform(sf, st_crs(to_crs))
  sv = vect(sf)
  ext(sv)
}

#' Prepare prediction table. Will predict values for every overpass in the
#' specified region and dates.
#'
#' @param features character vector of column names to train on
#' @param buffers_km vector of buffer radii in kilometers
#' @param hdf_root path to MCD19A2 HDF files
#' @param vrt_path directory to store overpass VRTs that reference HDF files.
#' @param this_date single row of pred_files data frame including the prediction
#'   date, year, and the corresponding XGBoost model file to use for predicting
#'   that date.
#' @param sat input "terra" or "aqua"
#' @param agg_level how much to aggregate the input MODIS AOD to use for large
#'   focal radii
#' @param agg_thresh use the aggregated AOD when `(radius * 1000 / agg_level /
#'   mcd_res) > agg_thresh`, where `radius` is in km, and `mcd_res` (resolution
#'   of input AOD) is in meters.
#' @param aoi sf shape of the region the model was trained over.
#' @param pred_bbox SpatExtent of the region to predict. If NULL, will crop
#'   using the shape provided in `aoi`.
pred_inputs <- function(features, buffers_km, hdf_root, vrt_path, load_sat,
                        this_date, agg_level, agg_thresh, aoi, pred_bbox = NULL){

  if('data.frame' %in% class(this_date)) this_date = this_date$pred_date

  # Get overpasses for the date
  if(length(this_date)>1) stop('Unexpected this_dates vector longer than 1')

  hdf_paths = list.files(file.path(hdf_root, format(this_date, '%Y.%m.%d')),
                         pattern = '\\.hdf$', full.names = TRUE)

  if(length(hdf_paths) == 0){
    stop('Unhandled case where all HDFs missing for date ', this_date)
  }

  binDT = bin_overpasses(hdf_paths)
  binDT = binDT[sat == load_sat]
  day_op = get_overpasses_vrts(hdf_paths, binDT, load_sat, vrt_path)

  # lapply by overpass, which are the groups in day_op
  tryCatch({
    day_rasters = mapply(FUN = rasterize_vrts,
                         op_vrts = day_op,
                         opDT = lapply(1:length(day_op), function(i) binDT[overpass_bin == i]),
                         SIMPLIFY = FALSE)
  }, error = function(e) {
    emsg = paste0('this_date = ', this_date, '\n',
                  'length(day_op) = ', length(day_op), '\n',
                  'binDT[, uniqueN(overpass_bin)] = ', binDT[, uniqueN(overpass_bin)], '\n',
                  'error = ', e)
    stop(emsg)
  })
  rm(day_op)

  aoi = sf_to_ext(aoi, crs_sinu)
  if(is.null(pred_bbox)) pred_bbox = aoi

  # Prepare predictors for each overpass
  op_to_table <- function(op_id){
    raslist = day_rasters[[op_id]]
    # note: hardcoded list item names (based on MCD19A2 resolution) and harcoded
    #   disaggregation factors
    r_relaz = terra::disagg(raslist$r4633, 5)
    r_times = terra::disagg(raslist$time, 1200)
    rstack = rast(list(raslist$r926, r_relaz, r_times))
    rm(r_relaz, r_times)
    # expand bbox to accomodate largest focal filter width
    # SpatExtent automatically subtracts the value on the minimum sides of the box, unlike st_bbox
    exp_bbox = pred_bbox + max(buffers_km)*1000
    rstack = terra::crop(rstack, exp_bbox)
    # keep a reference to the AOD layer for use in focal stats
    ras_aod = rstack$Optical_Depth_047
    aod_res = res(ras_aod)[1] # assuming square pixels
    # non-NA AOD layer
    ras_nna = !is.na(ras_aod)
    # aggregate AOD and non-NA for focal use
    agg_aod = terra::aggregate(ras_aod, fact = agg_level, fun = 'mean', na.rm = TRUE)
    agg_nna = terra::aggregate(ras_nna, fact = agg_level, fun = 'mean', na.rm = TRUE)
    # focal statistics
    focal_list = list()
    for(radius in buffers_km){
      if(radius * 1000 / agg_level / aod_res < agg_thresh){
        # use original resolution
        this_mat = circle_mat(radius, cellsize = aod_res, matrix = TRUE)
        foc_aod <- terra::focal(ras_aod, this_mat, fun = 'mean', na.rm = TRUE)
        foc_nna <- terra::focal(ras_nna, this_mat, fun = 'mean', na.rm = TRUE)
        focal_list[[paste0('Mean_AOD', radius, 'km')]]  <- foc_aod
        focal_list[[paste0('pNonNAAOD', radius, 'km')]] <- foc_nna
      } else {
        # use aggregated resolution
        this_mat = circle_mat(radius, cellsize = aod_res * agg_level, matrix = TRUE)
        foc_aod_agg <- terra::focal(agg_aod, this_mat, fun = 'mean', na.rm = TRUE)
        foc_nna_agg <- terra::focal(agg_nna, this_mat, fun = 'mean', na.rm = TRUE)
        fog_aod = terra::disagg(foc_aod_agg, agg_level)
        fog_nna = terra::disagg(foc_nna_agg, agg_level)
        rm(foc_nna_agg, foc_aod_agg)
        focal_list[[paste0('Mean_AOD', radius, 'km')]]  <- foc_aod
        focal_list[[paste0('pNonNAAOD', radius, 'km')]] <- foc_nna
      }
    }
    rm(radius)
    # stack focal results with the expanded crop above
    rstack = rast(list(rstack, rast(focal_list)))

    # crop to actual AOI
    rstack = terra::crop(rstack, pred_bbox)

    # convert to data.table
    rasDT = setDT(as.data.frame(rstack, na.rm = FALSE, xy = TRUE))
    rasDT = rasDT[!is.na(Optical_Depth_047)]
    setnames(rasDT, 'Optical_Depth_047', 'MCD19_AOD_470nm')
    rasDT[, overpass_time := as.POSIXct(overpass_time, tz = 'UTC', origin = '1970-01-01')]
    rasDT[, dayint := as.integer(as.Date(overpass_time))]

    # calc diffAOD for each buffer
    for(radius in buffers_km){
      rasDT[, c(paste0('diff_AOD', radius, 'km')) := MCD19_AOD_470nm - get(paste0('Mean_AOD', radius, 'km'))]
    }

    # convert AOD_QA to qa_best
    create_qc_vars(rasDT)

    # avoid error of writing NULL DT to FST
    if(nrow(rasDT) == 0) rasDT = data.table(NA)

    # overpass bin
    rasDT[, op_id := op_id]
  }
  all_ops = rbindlist(lapply(1:length(day_rasters), op_to_table), fill = TRUE)
  if('V1' %in% names(all_ops)) all_ops[, V1 := NULL]
  all_ops
}

#' Train a full model using DART and 2-fold CV to select hyperparameters from
#' maximin Latin hypercube sample
#'
dart_full <- function(
  data_train,
  y_var,
  features,
  n_rounds = 100,
  progress = TRUE
){
  xgb_threads <- get.threads()
  first_date = min(data_train$aer_date)
  data_train = data_train[, c(features, y_var), with = F]

  xdc_out <- xgboost.dart.cvtune(
    # by default, gives 100 rounds, and it is enough by experience
    n.rounds = n_rounds,
    d = data_train,
    dv = y_var,
    ivs = features,
    progress = progress,
    nthread = xgb_threads)

  # manually select and store some params
  param_dart <- c(xdc_out$model$params[c(1,2,4, 6:11)], nrounds = xdc_out$model$niter)

  # save the model
  # hash important input and output to create a unique name
  model_out_path = file.path('Intermediate',
                     paste0('full_model_dart_',
                            targets:::digest_obj64(list(data_train, features, y_var, param_dart)),
                            '.xgb'))
  xgboost::xgb.save(xdc_out$model, model_out_path)

  # record the prediction
  preds = xdc_out$pred.fun(data_train)
  #rmse_full <- sqrt(mean((preds - data_train[[y_var]])^2))

  # SHAP
  shap_pred <- as.data.table(xdc_out$pred.fun(data_train, predcontrib = TRUE, approxcontrib = FALSE))
  shap_bias <- first(shap_pred$BIAS)
  shap_pred[, BIAS := NULL]

  # features ranked by SHAP
  mean_shaps = colMeans(abs(shap_pred))
  ## feature names only:
  # features_rank_full_model <- names(mean_shaps)[order(mean_shaps, decreasing = TRUE)]
  # named vector:
  features_rank_full_model = mean_shaps[order(mean_shaps, decreasing = T)]

  list(features_rank_full_model = features_rank_full_model,
       y_preds = preds,
       shap_pred = shap_pred,
       shap_bias = shap_bias,
       model_out_path = model_out_path,
       first_date = first_date)
}

model_file_table = function(years, model_info_list){
  data.table(year = years, lapply(model_info_list, `[[`, 'model_out_path'))
}

#' Predict the difference between MCD19 and AERONET
#'
#' @param file_table is a data.frame containing the years, path to the model file, and prediction dates
run_preds = function(data, file_table, features){
  data = data[!is.na(MCD19_AOD_470nm)]
  data[, pred_date := as.Date(dayint, '1970-01-01')]

  # if the input data contains more than one model, run one model at a time
  outlist = list()
  for(pred_model in file_table[, unique(file_path)]){
    model = xgboost::xgb.load(pred_model)
    to_pred = data[pred_date %in% file_table[file_path == pred_model, pred_date], ]
    predvec = predict(model, as.matrix(to_pred[, features, with = FALSE]))
    # join predictions
    outpred = data.table(to_pred[, .(x, y, pred_date, op_id, MCD19_AOD_470nm)], preds = predvec)
    outpred[, MCD19_adjust := MCD19_AOD_470nm - preds]
    rm(to_pred)
    outlist[[length(outlist) + 1]] <- outpred
  }
  rbindlist(outlist, fill = TRUE)
}

# Maps ####

#' Compare adjusted AOD to original
#' @param data the data.table output of prediction
#' @param viz_date the single date to visualize from the output predictions
#' @param op_id the single overpass to visualize from the selected date
#' @return vertically stacked ggplots comparing original and adjusted MCD19A2 AOD
ggplot_orig_vs_adj = function(data, viz_date, viz_op){

  data = data[pred_date == viz_date & op_id == viz_op,
              .(x, y, MCD19_AOD_470nm, MCD19_adjust)]
  orig = simple.pred.map(data, fillvar = 'MCD19_AOD_470nm')
  adj = simple.pred.map(data, fillvar = 'MCD19_adjust')
  title = ggdraw() + draw_label(paste(as.character(viz_date), 'Overpass', viz_op))
  cowplot::plot_grid(title, orig, adj, ncol = 1, align = 'v')
}

# adapted from CONUS_air:plots.R
simple.pred.map = function(preds, fillvar, xvar = 'x', yvar = 'y',
                           scale.args = list()){
  ggplot(preds) +
  geom_raster(aes_string(xvar, yvar, fill = fillvar)) +
  do.call(scale_fill_distiller,
          c(list(palette = "Spectral"), scale.args)) +
  coord_equal() +
  theme_void() +
  theme(legend.justification = "left")
}

#' Interactively compare adjusted AOD to original
#'
#' @param data the data.table output of prediction
#' @param viz_date the single date to visualize from the output predictions
#' @param op_id the single overpass to visualize from the selected date
#' @param maxpixels resample raster version of predictions to approxmimately
#'   this many pixels and force the display in mapshot of this resolution
#' @return mapview object with a layer for the original and adjusted MCD19A2 AOD
mapshot_orig_vs_adj = function(data, viz_date, viz_op,
                               use_jenks = FALSE, maxpixels = NULL){

  data = data[pred_date == viz_date & op_id == viz_op,
              .(x, y, MCD19_AOD_470nm, MCD19_adjust)]
  ras = rasterFromXYZ(data, crs = crs_sinu)

  # resample raster
  if(!is.null(maxpixels)){
    mapviewOptions(georaster = TRUE, mapview.maxpixels = maxpixels)
    ras = sampleRegular(ras, maxpixels, asRaster = TRUE)
  }

  ub = unified_breaks(ras, n_classes = 10, viridis::inferno, use_jenks = use_jenks)
  get_mapview = function(i){
    mapview::mapview(ras[[i]], na.color = '#AAAAAA00', alpha.regions = 1,
                     col.regions = ub[[i]]$colors,
                     at = ub[[i]]$breaks, layer.name = names(ras[[i]]))}
  maps = lapply(1:2, get_mapview)
  if(!dir.exists(here('Intermediate/mapshot'))) dir.create(here('Intermediate/mapshot'), recursive = T)
  mapshot_path = here('Intermediate/mapshot',
                      paste0('orig_adj_MCD19_',
                             targets:::digest_obj64(list(data, use_jenks, maxpixels)),
                             '.html'))
  mapshot(maps[[1]] + maps[[2]], url = mapshot_path)
  out_size_MB = file.size(mapshot_path)/1024^2
  if(out_size_MB>500) warning('The output file', mapshot_path, 'is', round(out_size_MB,0), 'MiB in size')
  mapshot_path
}

#' Get unified breaks and color scheme for layer ranges that partially overlap
#' @param ras RasterStack or RasterBrick to prepare a unified set of class
#'   breaks for all its layers
#' @param n_classes number of classes, a numeric of length 1
#' @param color_func function to use to calculate colors, should take a single
#'   number as input
#' @return list with length equal to count of raster layers. Each item contains
#'   \code{$breaks}: numeric vector of class breaks, and \code{$colors} a
#'   character vector of hex-encoded colors.
unified_breaks = function(ras, n_classes, color_func, use_jenks = FALSE){
  rrange = range(ras[], na.rm = TRUE)
  if(use_jenks == TRUE){
    if(n_classes < 2) stop('cannot use n_classes < 2 with use_jenks = TRUE')
    # get all raster data into a single vector, removing NA values
    datavec = base::Reduce(c, as.data.table(ras[])[!is.na(get(names(ras)[1]))])
    nb1 = rgeoda::natural_breaks(k = n_classes, df = as.data.table(datavec))
    full_breaks = c(min(datavec), nb1, max(datavec))
    rm(datavec)
  } else {
    full_breaks = seq(rrange[1], rrange[2],
                      (rrange[2]-rrange[1])/n_classes)
  }
  colors = color_func(n_classes)

  breaks_by_layer = function(layernum, ras, full_breaks, colors){
    layer_min = min(ras[[layernum]][], na.rm = TRUE)
    layer_max = max(ras[[layernum]][], na.rm = TRUE)
    lyrI = intersect(which(full_breaks >= layer_min), which(full_breaks <= layer_max))
    lyrB = full_breaks[lyrI]

    if(full_breaks[min(lyrI)] > layer_min){
      lyrB = c(layer_min, lyrB)
      lyrI = c(min(lyrI)-1, lyrI)
    }
    if(full_breaks[max(lyrI)] < layer_max){
      lyrB = c(lyrB, layer_max)
      lyrI = c(lyrI, max(lyrI) + 1)
    }
    list(breaks = lyrB, colors = colors[lyrI[1:(length(lyrI)-1)]])
  }

  lapply(1:nlayers(ras), FUN = breaks_by_layer,
         ras = ras, full_breaks = full_breaks, colors)
}
