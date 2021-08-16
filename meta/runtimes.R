# See runtimes of previous runs
library(drake)
library(magrittr)
library(dplyr)
packageVersion("drake") # 7.8.0, Coco, before updating to R 4.0

allan_cache = drake_cache(path = "/scratch/cache/aeronet_drake")      # Dec 22 2019, Allan's run
allan_cache = drake_cache(path = "/data-coco/mcd19/cache/aeronet_drake") # Jan 4 2020, newer run
# changed group to shared and added group write to everything; deleted drake/history/lock
yang_cache  = drake_cache(path = "/scratch/cache/aeronet_drake_yang") # Aug 28 2019
# Johnathan's old cache on Belle /scratch/cache/aeronet_drake, but only has 13 targets, and was probably just aeronet_subset project that Yang built on
john_cache = drake_cache(path = "/data-belle/mcd19/cache/aeronet_drake_johnathan")

# note: should add years to target names (like in MAIAC processing) unless we do not want to keep cache of previous years

yc = cached(cache = yang_cache)
ac = cached(cache = allan_cache)
jc = cached(cache = john_cache)
length(ac) # 2956
length(yc) # 2955
length(jc) # 14

abt = build_times(cache = allan_cache, all_of(ac))
ybt = build_times(cache = yang_cache, all_of(yc))
jbt = build_times(cache = john_cache, all_of(jc))

# checking time to run aod_MODIS_newVars() for Allan's run:
library(dplyr)
abt %>% filter(substr(target, 1, 7) == 'new_var') %>% select(elapsed) %>% summary()
#     elapsed                    
# Min.   :0.268s                
# 1st Qu.:34.2535s              
# Median :38.215s               
# Mean   :46.5202421340629s     
# 3rd Qu.:43.8105s              
# Max.   :119s (~1.98 minutes) 
abt %>% filter(substr(target, 1, 7) == 'new_var') %>% select(elapsed) %>% sum()
# 34006 seconds, or 9.4 hours for CONUS

library(fst)
write_fst(abt, "meta/allan_runtimes.fst", compress = 100)
write_fst(ybt, "meta/yang_runtimes.fst", compress = 100)
#write_fst(jbt, "meta/johnathan_runtimes_pre-optimize.fst", compress = 100)
#write_fst(jbt, "meta/johnathan_runtimes.fst", compress = 100)

library(dplyr)
library(stringr)
abt %>% arrange(desc(elapsed))
ybt %>% arrange(desc(elapsed))

sum(abt$elapsed)/3600

# longest is ref_sel_aer at 2.1 or 1.7 hours. 
# wPred and aer_data are around 10 min. 
# next is mcd_DAY_aqua/terra: each of these around 4.5 minutes for Allan or 9 minutes for Yang
abt %>% select(target, elapsed) %>% arrange(desc(elapsed)) %>% head(20)
ybt %>% select(target, elapsed) %>% arrange(desc(elapsed)) %>% head(20)
# 2x difference betwen Allan and Yang runtimes. Allan likely used used 1/2 as many simultaneous processes, Yang may have used all cores including HT cores
# yes: user time is 2x for Allan, so data.table had access to more threads during his run
# but the system overhead was also 2x, so overall runtime was probably longer for Allan
abt %>% filter(str_detect(target, pattern = "^mcd_")) %>% arrange(desc(elapsed))
ybt %>% filter(str_detect(target, pattern = "^mcd_")) %>% arrange(desc(elapsed))

# cv_results_aqua/terra are fast because they parallelize across all cores
# aer_data has 2.5x user: could it be given more threads? Has 1/2 system vs. elapsed, so may be IO limited. 
abt %>% arrange(desc(user))
ybt %>% arrange(desc(user))

# mostly just see the 2x system time for Allan's run of mcd_DAY_aqua/terra
abt %>% arrange(desc(system))
ybt %>% arrange(desc(system))

