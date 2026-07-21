#!/bin/bash

# Predict corrected MCD19A2 AOD values for every cell and day.

#SBATCH --job-name Earth_obs_cleaning
#SBATCH --output eoc_%04a.txt
#SBATCH --partition batch
#SBATCH --time 4:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 64G
#SBATCH --array 1-1240

D=/oscar/data/allanjust/Earth_obs_cleaning
cd $D/src
apptainer exec \
    -B $D:/data \
    $D/Earth_obs_cleaning.simg \
    R -q --slave -e \
    "source('code/slurm.R'); make.preds($SLURM_ARRAY_TASK_ID)"
echo 'Done'
