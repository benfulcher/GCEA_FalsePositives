#!/bin/bash

# First compute null distributions under null maps:

# Spatially autocorrelated:
sbatch slurm_mouseCategoryNullCompute_spatialLag_wholeBrain.sh
sbatch slurm_mouseCategoryNullCompute_spatialLag_cortex.sh

# Independent random:
sbatch slurm_mouseCategoryNullCompute_randomMap_wholeBrain.sh
sbatch slurm_mouseCategoryNullCompute_randomMap_cortex.sh

# Intra-category correlation nulls (not really used anymore):
# sbatch slurm_mouseIntraCorr_raw.sh

# Compute all of my false-positive significance rates (under conventional enrichment methodology):
sbatch slurm_mouseUniformRandomShuffleUniformRandom.sh
sbatch slurm_mouseSpatialLag.sh
sbatch slurm_mouseUniformRandom.sh
