# import AERONET data
library(drake)
library(here)
library(data.table)
library(purrr)
library(sp)
library(maptools)
library(dplyr)
library(ggplot2)
library(ggmap)
library(raster)

# https://aeronet.gsfc.nasa.gov/aeronet_locations_v3.txt

# path to downloaded Aeronet data (e.g. version 3 level 2.0)
# https://aeronet.gsfc.nasa.gov/cgi-bin/combined_data_access_v3

plan <- drake_plan(
  aeronet_stations = fread(file = file_in(here("data/aeronet_locations_v3.txt"))),
  setnames(aeronet_stations, c("Site_Name", "lon", "lat", "elevm")),
  
  aeronet = fread(file_in('unzip -cq data/All_Sites_Times_All_Points_AOD20.zip')),
  aeronet[, day := as.Date(strptime(`Date(dd:mm:yyyy)`, format = "%d:%m:%Y"))],
  # pryr::object_size(aeronet)
  # nrow(aeronet)
  # subset by date
  daterange = data.table(start = as.Date("2015-01-01"), end = as.Date("2018-11-01")),
  aeronet_sub = aeronet[day >= daterange$start & day <= daterange$end, ],
  # pryr::object_size(aeronet)
  nrow(aeronet_sub)
  
)
plan
make(plan)

config <- drake_config(plan)
vis_drake_graph(config)



aeronet[, day := as.Date(strptime(`Date(dd:mm:yyyy)`, format = "%d:%m:%Y"))]
# pryr::object_size(aeronet)
# nrow(aeronet)
# subset by date
daterange = data.table(start = as.Date("2015-01-01"), end = as.Date("2018-11-01"))
aeronet = aeronet[day >= daterange$start & day <= daterange$end, ]
# pryr::object_size(aeronet)
# nrow(aeronet)


setnames(aeronet, "Precipitable_Water(cm)", "aeronet_cwv") 


##### old code ####

aeronetfiles

test <- fread("data/aeronet/AOT/LEV20/ALL_POINTS/920801_170325_GSFC.lev20")
test
dput(names(test))
gsfc <- test[, c("Date(dd-mm-yy)", "Time(hh:mm:ss)", "Julian_Day", "AOT_440", "Water(cm)", "Solar_Zenith_Angle"), with = F]
gsfc[, day := as.Date(strptime(`Date(dd-mm-yy)`, format = "%d:%m:%Y"))]
summary(gsfc$day)
# delete old date field
gsfc[, `Date(dd-mm-yy)` := NULL]


table(aeronetfiles$Site_Name %in% aeronetsites$Site_Name)
table(aeronetsites$Site_Name %in% aeronetfiles$Site_Name)
# inner join (sites with data and coordinates) - there are some site names without data files
aerojoin <- merge(aeronetsites, aeronetfiles, by = "Site_Name")
aerojoin

# I want to restrict aerojoin to the points that fall in roi
coordinates(aerojoin) <- ~ lon + lat # makes aerojoin a SPDF
# let's load roi as a SpatialPolygonsDataFrame
roi <- readShapePoly("data/gis/nemia.shp") # warning to use readOGR or st_read
# let's plot 
plot(roi)
plot(aerojoin, add = TRUE)
# let's restrict to the points inside roi
nrow(aerojoin)  # 973
aerojoin <- aerojoin[!is.na(over(aerojoin, roi)$OBJECTID), ]
nrow(aerojoin)  #  90
# let's plot again
plot(roi)
plot(aerojoin, add = TRUE)
# let's go back to data.table
aerojoin <- as.data.table(aerojoin)

# let's also check with ggmap 
extentroi <- as.vector(extent(roi))
mymap <- get_map(c(extentroi[1], extentroi[3], extentroi[2], extentroi[4]), 
                 source = "stamen", maptype = "terrain-background")
ggmap(mymap) +
  geom_point(aes(lon, lat), aerojoin) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank())

# let's import all the files with AOT in roi
importfile <- paste0("data/aeronet/AOT/LEV20/ALL_POINTS/", as.character(aerojoin$file))
# I also restrict to the columns I want
aotdt <- importfile %>% 
  parallel::mclapply(fread) %>% 
  parallel::mclapply(function(x) x[, c("Date(dd-mm-yy)", "Time(hh:mm:ss)", "Julian_Day", "AOT_440", "Water(cm)", "Solar_Zenith_Angle"), with = FALSE]) 
#aotdt <- importfile %>% 
#  map(fread) %>% 
#  map(function(x) x[, c("Date(dd-mm-yy)", "Time(hh:mm:ss)", "Julian_Day", "AOT_440", "Water(cm)", "Solar_Zenith_Angle"), with = FALSE]) 
  # this is pretty slow. at least partially from j drive. possibly copy files to server, possibly replace DT step's map with mclapply or futures
  # lots of warnings about changing type to character when encountering 'N/A' values. Can set 'colClasses' to character on loading to avoid issues. 
# let's set names to the list and rbind them 
names(aotdt) <- aerojoin$Site_Name
aotdt <- rbindlist(aotdt, idcol = "Site_Name") 
str(aotdt)
aotdt[, day := as.Date(strptime(`Date(dd-mm-yy)`, format = "%d:%m:%Y"))]
aotdt[, `Date(dd-mm-yy)` := NULL] # delete old date field
# let's define mins from midnight 
aotdt[, mins := lapply(`Time(hh:mm:ss)`, function(x) {xnum <- as.numeric(strsplit(x, ":")[[1]]); xnum[1]*60 + xnum[2] + (xnum[3]/60)}), by = 1:nrow(aotdt)]
str(aotdt)
# I guess we want AOT_440, and CWV to be numeric
aotdt[, `:=`(AOT_440 = as.numeric(AOT_440), `Water(cm)` = as.numeric(`Water(cm)`))]
str(aotdt)
# are there missing values for data columns?
nrow(aotdt)                               #  795033
nrow(aotdt[is.na(AOT_440)])               #    1707
nrow(aotdt[is.na(`Water(cm)`)])           #   22541
nrow(aotdt[is.na(Solar_Zenith_Angle)])    #       0

# let's drop them 
aotdt <- aotdt[!is.na(AOT_440)] 
nrow(aotdt)                               #  793326
# keeping NAs for two new columns

# let's restrict to the time period we are interested in 
aotdt <- aotdt[day >= as.Date("2000-01-01")]
aotdt <- merge(aotdt, aerojoin[, .(Site_Name, lon, lat)], by = "Site_Name", all.x = TRUE)
nrow(aotdt)                               #  737915

# let's save this dataset in intermediate 
  saveRDS(aotdt, file = paste0("data/intermediate/AOT_CWV_SZE_aeronet_v", Sys.Date(), ".rds"))
# also FST - might make reading from J drive quicker 
  fst::write_fst(aotdt, paste0("data/intermediate/AOT_CWV_SZE_aeronet_v", Sys.Date(), ".fst"), compress = 90) # much faster, bit bigger than .rds

# how many stations per year? 
aotdt[, length(unique(Site_Name)), keyby = format(day, "%Y")]   # range of 6 to 51


# skipping below for CWV run
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # minimal dataset to send to Itai 
# # we care only after 2000
# aotshort <- aotdt[, list(minday = min(day),
#                          maxday = max(day)), by = c("Site_Name", "lon", "lat")]
# # we want to save a dt to send to Itai
# # write.csv(aotshort, file = paste0("data/aeronet/aeronet_stations_in_nemia_v", Sys.Date(), ".csv"), row.names = FALSE)
# 
# 
# # to send to Alex 
# library(data.table)
# aotshort <- fread("data/aeronet/aeronet_stations_in_nemia_v2017-03-29.csv")
# aotshort[, c("minday", "maxday") := lapply(.SD, as.Date), .SDcols = c("minday", "maxday")]
# aotlong <- aotshort[, list(day = seq.Date(from = minday, to = maxday, by = 1), 
#                            lon = lon,
#                            lat = lat), by = Site_Name]
# aotlong
# 
# # end file 