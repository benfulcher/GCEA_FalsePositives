#!/bin/bash

# ------------------------------------------------------
# Compute null distributions under different phenotype ensembles:
# ------------------------------------------------------
sbatch slurm_humanCategoryNullCompute_randomMap.sh
sbatch slurm_humanCategoryNullCompute_spatialLag.sh

# ------------------------------------------------------
# Compute false-positive significance rates (under conventional GO enrichment):
# ------------------------------------------------------
sbatch slurm_humanUniformRandomShuffleUniformRandom.sh
sbatch slurm_humanUniformRandom.sh
sbatch slurm_humanSpatialLag.sh

# Intra-category correlation (but randomized versions are not used):
sbatch slurm_humanIntraCorr_raw.sh
