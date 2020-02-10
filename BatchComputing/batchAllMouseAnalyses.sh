#!/bin/bash

# ------------------------------------------------------
# First compute null distributions under (null) phenotype ensembles:
# ------------------------------------------------------

# Spatially autocorrelated:
sbatch slurm_mouseCategoryNullCompute_spatialLag_wholeBrain.sh
sbatch slurm_mouseCategoryNullCompute_spatialLag_cortex.sh

# Independent random:
sbatch slurm_mouseCategoryNullCompute_randomMap_wholeBrain.sh
sbatch slurm_mouseCategoryNullCompute_randomMap_cortex.sh

# ------------------------------------------------------
# Compute FPSR (under conventional GSEA):
# ------------------------------------------------------
sbatch slurm_mouseUniformRandomShuffleUniformRandom.sh
sbatch slurm_mouseUniformRandom.sh
sbatch slurm_mouseSpatialLag.sh

# Intra-category correlation nulls (randomized versions are no longer used):
sbatch slurm_mouseIntraCorr_raw.sh
