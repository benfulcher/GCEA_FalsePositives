#!/bin/bash

# Compute null distributions under different phenotype ensembles:
sbatch slurm_humanCategoryNullCompute_randomMap.sh
sbatch slurm_humanCategoryNullCompute_spatialLag.sh

# Intra-category correlation nulls (not really used anymore):
# sbatch slurm_humanIntraCorr_raw.sh

# Compute all false-positive significance rates (under conventional enrichment methodology)
sbatch slurm_humanUniformRandomShuffleUniformRandom.sh
sbatch slurm_humanSpatialLag.sh
sbatch slurm_humanUniformRandom.sh
