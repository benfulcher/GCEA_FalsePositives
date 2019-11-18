#!/bin/bash

# First compute null distributions under different maps (20k null samples per category):
sbatch slurm_mouseCategoryNullCompute_spatialLag_wholeBrain.sh
sbatch slurm_mouseCategoryNullCompute_spatialLag_cortex.sh

sbatch slurm_mouseCategoryNullCompute_randomMap_wholeBrain.sh
sbatch slurm_mouseCategoryNullCompute_randomMap_cortex.sh

# Intra-category correlation nulls (not really used anymore):
sbatch slurm_mouseIntraCorr_VE1.sh
sbatch slurm_mouseIntraCorr_raw.sh

# Compute all of my false-positive proportions (across 10k random maps):
sbatch slurm_mouseUniformRandomShuffleUniformRandom.sh
sbatch slurm_mouseSpatialLag.sh
sbatch slurm_mouseUniformRandom.sh
