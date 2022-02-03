# This runs the workflow in _targets.R in stages, running some parts in
# parallel using future, and others parts sequentially.
# Portions running XGBoost are run sequentially because it uses OpenMP for
# parallelism, and running multiple XGBoost jobs would overload the available cores.

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
future_workers = 10L
end_part1 = c('traindata_aqua_conus', 'traindata_terra_conus')
message('Part 1: Preparing training data using ', future_workers, ' workers at ',
        Sys.time())
tar_make_future(names = !!end_part1, workers = future_workers)
check_progress(1, end_part1)

# 2. Train full model and Generate CV report for all years, serial ####
message('Part 2: Generating CV report and training full models at ', Sys.time())
end_part2 = c('initial_cv_report', 'full_model_aqua_conus', 'full_model_terra_conus')
tar_make(names = !!end_part2)
check_progress(2, end_part2)

# 3. Prepare CONUS prediction inputs, parallel ####
# Predictions run for dates in `pred_dates`.
# An error will be thrown if trying to predict for a date in a year that has not
# had a model trained.
message('Part 3: Preparing CONUS prediction inputs using ', future_workers,
        ' workers at ', Sys.time())
end_part3 = c('predinput_aqua_conus', 'predinput_terra_conus')
tar_make_future(names = !!end_part3, workers = future_workers)
check_progress(3, end_part3)

# 4. Run predictions, serial ####
message('Part 4: Running CONUS predictions at ', Sys.time())
end_part4 = c('pred_out_aqua_conus', 'pred_out_terra_conus')
tar_make(names = !!end_part4)
#check_progress(4, end_part4)
message('Finished part 4 at ', Sys.time())

# 5. Visulize predictions ####
