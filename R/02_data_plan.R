
data_plan <- drake_plan(
  # 1.  Aeronet ---------------------------------------------------------------------
  
  # AERONET station locations file: 
  aer_stns = fread(file = aer_stn_path, col.names = c("Site_Name", "lon", "lat", "elevm")),
  region_buff = target(get_aoi_buffer(aoiname), transform = map(aoiname = !!aoiname)),
  aerpts = get_aer_spatial(aer_stns),               # wgs84
  aer = target(select_points(aerpts, region_buff),
               transform = map(aoiname)),           # AERONET sites in region as points
  
  # get nearest MODIS cells 
  refgrid = get_ref_grid(refgrid_path),
  
  # aeronet sites with cloest MODIS grid
  nearest_grid = target(get_nearest_cell(aer, refgrid), transform = map(aer)), # Aeronet matched to grid
  nearest_grid_2 = target(remove_site_on_water(aer_nearest), transform = map(aer)), # remove 2 sites on water
  
  # limit aeronet data by date
  aer_data = target(get_stn_data(aod_dir = aer_files_path, stn_names = aer_nearest_2$Site_Name,
                                 date_start, date_end),  # 26171206x56, 11G
                    transform = map(date_start = !!date_start, 
                                    date_end = !!date_end),
                    format = "fst_dt"),
  # aer_btw = target(sel_data_bytime(aer_data, date_start, date_end),
  #                  transform = map(date_start = !!date_start, 
  #                                  date_end = !!date_end,
  #                                  .id = FALSE)),
  # # limit aeronet data by station, could be NEMIA or CONUS 
  # sel_aer_region = target(sel_data_bystation(aer_btw, aer_nearest_2),
  #                         transform = map(aer_btw, .id = FALSE)),
  
  
  # 2. Interpolation -----------------------------------------------------------
  # so instead of running all the Aeronet data, only run what is used
  # it also join with the nearest grid id 
  wPred = target(interpolate_aod(aer_data, nearest_grid_2), 
                 transform = map(aer_data, nearest_grid_2),
                 format = "fst_dt"),
  
  # 3. WPred Join MCD19 -----------------------------------------------------------
  # mcd is the large MCD19 dataset, limit to `refsub`:nearby(300km) when readin to 
  # control total lines, unable to read in all data together, too many rows
  mcd = target(read_mcd19_one(i, sat = sat, mcd19path_CONUS), 
               transform = cross( i = !!MCD_files_i, 
                                  sat = !!sats),
               format = "fst_dt"),
  # join is the rolling-joined dataset with aeronet
  join = target(rolling_join(mcd_sat = mcd, aer_data_wPred = wPred),
                transform = cross(mcd, wPred, .id = mcd),
                format = "fst_dt"),
  
  # For calculating new variables
  # selected aeronet 
  sel_aer = target(get_conus_aer_used(aer, unique(aer_data$AERONET_Site_Name)),
                   transform = map(aer, aer_data, .id = FALSE)),
  
  # buffers 
  ref = target(ref_in_buffer(sel_aer, refgrid),
               transform = map(sel_aer)),
  # sites
  p = target(get_site_list(sel_aer), transform = map(sel_aer)),
  n = target(get_site_name(sel_aer), transform = map(sel_aer)),
  
  # mapping 
  # this way doesn't work:
  # sub_grid = target(select_refgrid_subset(ref_list0 = ref, point0 = p, name0 = n),
  #              transform = map(ref, p, n))
  sub_grid = target(purrr_pmap(ref = ref, p = p, n = n),
                    transform = map(ref, p, n, .id = FALSE),
                    format = "fst_dt"),
  # calculate new vars to obtain the dataset ready for modeling for terra and aqua, using join, mcd, sub_grid 
  new_var = target(aod_MODIS_newVars(aod_join_MODIS = join, 
                                     MODIS_all = mcd, refDT_sub = sub_grid),
                   transform = map(join, mcd, .id = c(join)),
                   format = "fst_dt"),
  # deliver two rowbinded dataset: for terra and aqua 
  # I am not sure what this way doesn't work: 
  # binded = target(rbindlist(new_var),
  #                 transform = combine(new_var, .by = c(i, sat))),
  # `.id_chr` let you save the file using target names:
  out = target(write_new_var(new_var, id = .id_chr), transform = map(new_var))
)
