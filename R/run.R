# This will run the workflow in _targets.R in stages:
# The first part has parallelism managed by targets and future.
# The second part for models is run sequentially by targets, since XGBoost
# uses OpenMP and all cores by default to parallelize work.

# Note the use of `shortcut = TRUE` means we're relying on checks against
# tar_outdated() to avoid mismatched dependencies

library(targets)

# 1. Prepare predictors ####
end_part1 = 'modelinput_'
message('Starting part 1 using parallel futures at ', Sys.time())
tar_make_future(names = starts_with(!!end_part1),
                workers = 8L)
message('Finished part 1 at ', Sys.time())

# check if first part has completed successfully
message('Checking status of part 1 before starting part 2...')
outdated = tar_outdated()
unfinished = grep(end_part1, outdated)
if(any(unfinished)){
  stop('Aborting workflow. The following targets from part 1 need to complete before running part 2:\n',
       paste('  ', sort(outdated[unfinished]), collapse = '\n'))}

# 2. Train models ####
message('Starting part 2 sequentially at ', Sys.time())
end_part2 = 'initial_cv_'
tar_make(names = starts_with(!!end_part2), shortcut = TRUE)
message('Finished part 2 at ', Sys.time())

# check if training was successful
message('Checking status of part 2 before starting part 3...')
outdated = tar_outdated()
unfinished = grep(end_part2, outdated)
if(any(unfinished)){
  stop('Aborting workflow. The following targets from part 2 need to complete before running part 3:\n',
       paste('  ', sort(outdated[unfinished]), collapse = '\n'))}

# 3. Render report ####
message('Rendering report at ', Sys.time())
tar_make(names = initial_cv_report, shortcut = TRUE)
message('Done')
