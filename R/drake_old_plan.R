nemia_plan <- drake_plan(
  
  
  # AERONET station locations file: 
  aer_stns = fread(file = file_in(aer_stn_path), col.names = c("Site_Name", "lon", "lat", "elevm")),
  conus_buff = get_conus_buff(),
  nemia_buff = get_nemia_buff(),
  aerpts = get_aer_spatial(aer_stns),               # wgs84
  conus_aer = select_points(aerpts, conus_buff),    # AERONET sites in CONUS as points
  nemia_aer = select_points(conus_aer, nemia_buff), # AERONET sites in NEMIA as points
  
  # get nearest MODIS cells 
  refgrid = get_ref_grid(),
  candidate_refpts = points_in_buffer(nemia_aer, refgrid),
  aer_nearest = get_nearest_cell(nemia_aer, candidate_refpts), # Aeronet matched to grid, 174x4
  aer_nearest_NEMIA = remove_site_on_water(aer_nearest), # remove 2 sites 
  
  # measurements. Editing values of date_start and date_end should make use of the cache
  # if range was previously calc'd
  aer_data = get_stn_data(file_in(aer_data_path), aer_columns),  # 26171178x6
  aer_btw = target(sel_data_bytime(aer_data, date_start, date_end),
                   transform = map(date_start = !!date_start, date_end = !!date_end)),
  
  # selected measurements, transform by same date range if not automatic in name
  nemia_data = target(sel_data_bystation(aer_btw, nemia_aer),
                      transform = map(aer_btw)), # use map to map to aer_btw
  
  # joined mcd19 MAIAC values. transform by same date range if not automatic in name
  
)
