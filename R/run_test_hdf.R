# Test HDF loading in training side. Memory demands are high for CONUS, perhaps
# up to 90GB per worker.
# This runs the workflow in _targets.R in stages, running some parts
# in parallel using future, and others parts sequentially.

library(targets)
Sys.setenv(RSTUDIO_PANDOC = '/usr/lib/rstudio-server/bin/pandoc')

# Check the progress of the current batch of targets before starting the next.
check_progress = function(finished_part, check_targets){
  message('Finished part ', finished_part, ' at ', Sys.time())
  message('Checking status of part ', finished_part, ' before starting part ',
          finished_part + 1, '...')
  outdated = tar_outdated(!!check_targets)

  if(length(outdated) > 0){
    stop('The following targets from part ', finished_part,
         ' need to complete before running part ', finished_part + 1, ':\n',
         paste('  ', sort(outdated), collapse = '\n'))
  } else {
    message('Part ', finished_part, ' complete.')
  }
}

# 1. Prepare training data, parallel ####
# training runs for all years in `process_years`.
future_workers = 5L # current HDF loading and merging has high memory demands
end_part1 = c('traindata_aqua_conus')
#end_part1 = c('traindata_aqua_conus', 'traindata_terra_conus')
message('Part 1: Preparing training data using ', future_workers, ' workers at ',
        Sys.time())
tar_make_future(names = !!end_part1, workers = future_workers)
check_progress(1, end_part1)

# # 2. Generate CV report for all years, serial ####
# message('Part 2: Generating CV report at ', Sys.time())
# end_part2 = c('initial_cv_report')
# tar_make(names = !!end_part2)
# check_progress(2, end_part2)
