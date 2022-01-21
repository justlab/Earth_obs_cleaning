# This will run the workflow in _targets.R in stages:
# The first part has parallelism managed by targets and future.

# The second part for models is run sequentially by targets, since XGBoost
# uses OpenMP and all cores by default to parallelize work.

library(targets)

# 1. Prepare predictors ####
# mcd19_vars_aqua_conus prepares IVs for the training side, runs for all years in `process_years`
# predinput_aqua_conus prepares IVs for the prediction side, runs for dates in `pred_dates`
end_part1 = c('mcd19_vars_aqua_conus', 'predinput_aqua_conus')
message('Starting part 1 using parallel futures at ', Sys.time())

tar_make_future(names = !!end_part1,
                workers = 8L)
message('Finished part 1 at ', Sys.time())

# check if first part has completed successfully
message('Checking status of part 1 before starting part 2...')
outdated = tar_outdated(!!end_part1)
unfinished = grep(paste0(end_part1, collapse = '|'), outdated)

if(any(unfinished)){
  stop('Aborting workflow. The following targets from part 1 need to complete before running part 2:\n',
       paste('  ', sort(outdated[unfinished]), collapse = '\n'))}

# 2. Train full model, prepare CONUS prediction table, run predictions ####
message('Starting part 2 sequentially at ', Sys.time())
tar_make(names = pred_out_aqua_conus)
message('Finished part 2 at ', Sys.time())

# # 3. Generate CV report for all years
# tar_make(names = initial_cv_report)
