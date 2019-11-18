#!/bin/bash

# First compute null distributions under different maps (20k null samples per category):
sbatch slurm_humanCategoryNullCompute_randomMap.sh
sbatch slurm_humanCategoryNullCompute_spatialLag.sh

# Intra-category correlation nulls (not really used anymore):
sbatch slurm_humanIntraCorr_VE1.sh
sbatch slurm_humanIntraCorr_raw.sh

# Compute all false-positive significance rates:
sbatch slurm_humanUniformRandomShuffleUniformRandom.sh
sbatch slurm_humanSpatialLag.sh
sbatch slurm_humanUniformRandom.sh
