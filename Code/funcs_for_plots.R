# Functions for AOD
# Yang Liu


# for data ----------------------------------------------------------------


read.lst <- function(sat = "terra", columns0 = NULL){
  if (sat != "terra") choose = "mcd19_T" else choose = "mcd19_A" # load terra by default
  lst_files <- list.files(path = "/data-coco/mcd19/fst/nemia/2018", pattern = choose, full.names = T)
  dt = rbindlist(lapply(lst_files, function(x)read.fst(x, as.data.table = T, columns  = columns0)))
  if(is.null(columns0)) dt[,c("x", "y", "inNEMIA"):=NULL]
  
  dt[, join_time:=overpass_time]
  setkey(dt, idLSTpair0, join_time)
  return(dt)
}



#  general functions for making plots -------------------------------------


# Scatterplot -------------------------------------------------------------

plot.scatter <-  function(data, x, y, size0 = 0.1, alpha0 = 0.3, dilute = F, dilute_by = 10){
  set.seed(1234)
  if (dilute) {
    data <- data[sample(nrow(data), nrow(data)/dilute_by)] # dilute
  }
  plot0 <- ggplot(data = data, aes(x = data[[x]], y = data[[y]]))+
    # geom_point(size = size0, alpha = alpha0) + 
    geom_density_2d(aes(fill = ..level..), geom = "polygon") + 
    labs(x = x, y = y) + 
    theme_bw()
  # add histogram
  ggExtra::ggMarginal(plot0, type = "histogram", bins = 50, size = 10, color="white")
}

 
# Maps --------------------------------------------------------------------
# us map using ggmap
plot.us.ggmap <- function(data0){
  # data0 contains the points to plot, plot on whole us mainland
  map_wholeUS <- readRDS("~/Jdrive/PM/Just_Lab/projects/ECHO_PM/data/intermediate/map_us_all.RData")
  ggmap(map_wholeUS) + 
    geom_point(data = data0, aes(x = lon, y = lat),size = 1, color= "red") + # MOD
    # geom_label(data = SiteName_aodid, aes(label = Site_Name), size = 2, alpha = 0.7) + 
    theme(legend.position='none') 
}


# us map with fill
plot.us.map.size <- function(data0 = aeronetdt, size0 = "Ndays", sat0 = sat, NEMIA = F){
  # optional to print NEMIA alone
  NE_list <- c("maine", "new hampshire", "vermont",
               "massachusetts", "rhode island", "connecticut", 
               "new york", "new jersey", "pennsylvania",
               "delaware", "maryland", "west virginia", "virginia")
  
  # borrowed from nemia_aeronet_stn.R
  g1 = ggplot(data = if (NEMIA) map_data("state", region = NE_list) else map_data("state"), 
              aes(x = long, y = lat, group = group)) +
    geom_polygon(fill = 'NA', color = "grey") + 
    coord_map("albers", parameters = list(at0 = 45.5, lat1 = 29.5), expand = c(0, 0)) +
    geom_point(aes_string("lon", "lat", "group" = 1, size = size0),
               alpha = 0.3, data = data0)+
    # geom_point(aes(lon, lat, group = 1, size = stn_count), 
    #            alpha = 0.3, data = data0)+
    scale_size_continuous(name = paste0("Number of \nDays (",Hmisc::capitalize(sat0),")"),
                          breaks = c((1:3)*500), range = c(1, 10)) +
    ylab("Latitude") + xlab("Longitude") + 
    theme_minimal() +
    theme(legend.position = c(0.1, 0.75))
  if (!NEMIA){
    g1 <- g1 +
      xlim(c(-130, -60)) +
      ylim(c(25, 50))
  }
  g1
}


plot.nemia.map <- function(data0){
  # plot NEMIA region (NY state)
  map_us <- readRDS("~/Jdrive/PM/Just_Lab/projects/ECHO_PM/data/intermediate/map_us_NE.RData")
  ggmap(map_us) + 
    geom_point(data = data0, aes(x = lon, y = lat),size = 0.05) + # MOD
    # geom_label(data = SiteName_aodid, aes(label = Site_Name), size = 2, alpha = 0.7) + 
    theme(legend.position='none') 
}

plot.nyc <- function(data, lon = 'lon', lat = 'lat', zoom_in_option = 2, fill = NULL){
  # only need lat lon (WGS84) from data, will plot all the points 
  
  library(data.table); library(fst); library(sf); library(ggplot2)
  nemia_poly <- readRDS("~/Jdrive/PM/Just_Lab/projects/ECHO_PM/Intermediate/nemia_poly_WGS.rds")
  data <- as.data.table(data)
  setnames(data, c(lon, lat), c("lon", "lat"))
  
  bd1 <- as.data.frame(cbind(
    lon = c(-74.4, -73.5),
    lat = c(41.5, 40.5)
  ))
  
  # zoom in 
  bd2 <- as.data.frame(cbind(
    lon = c(-74.25, -73.7),
    lat = c(40.98, 40.4)
  ))
  # zoom out 
  bd3 <- as.data.frame(cbind(
    lon = c(-75, -73),
    lat = c(42, 40)
  ))
  bd <- as.matrix(cbind(
    ny1 = c((bd1)[1, 1], (bd1)[2, 1], (bd1)[2, 2], (bd1)[1, 2]),
    ny2 = c((bd2)[1, 1], (bd2)[2, 1], (bd2)[2, 2], (bd2)[1, 2]),
    ny3 = c((bd3)[1, 1], (bd3)[2, 1], (bd3)[2, 2], (bd3)[1, 2])
  ))
  
  fi = zoom_in_option
  data_nyc <- data[lon>bd[1,fi] & lon<bd[2,fi] & lat>bd[3,fi] & lat< (bd[4,fi]), ]
  nycmap <- ggplot() + 
    # NEMIA land background
    geom_sf(data = nemia_poly, fill = "gray90", color="black", inherit.aes = F) +
    # map extent
    xlim(bd[1,fi], bd[2,fi]) +
    ylim(bd[3,fi], bd[4,fi]) + 
    # optional: overlay NEMIA outline on top of cells
    geom_sf(data = nemia_poly, fill = NA, color="black", inherit.aes = F) +
    # # scale use ggsn
    # ggsn::scalebar(dist = 6, transform = FALSE, dist_unit = 'km', 
    #         border.size = 0.4, st.size = 4.5, # boarder width and font size
    #         anchor = c(x = 1761000, y = 2224000),
    #         box.fill = c('black','white'),
    #         x.min = bd[1,fi], x.max = bd[2,fi], y.min= bd[3,fi], y.max =  bd[4,fi]) +
    coord_sf(expand = F) + 
    theme_bw() 
  if(is.null(fill)){
    nycmap <-  nycmap + geom_point(data = data_nyc, aes(x = lon, y = lat))
  } else {
    nycmap <-  
      nycmap + geom_point (data = data_nyc, aes_string(x = "lon", y = "lat", color = fill))+
      scale_color_viridis_c(option = "plasma", direction = -1)
    
  }
  return(nycmap)
}
