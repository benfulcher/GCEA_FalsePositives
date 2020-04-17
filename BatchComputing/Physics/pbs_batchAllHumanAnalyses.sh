#!/bin/bash

# ------------------------------------------------------
# Compute null distributions under different phenotype ensembles:
# ------------------------------------------------------
qsub pbs_humanCategoryNullCompute_SBP-rand.sh
qsub pbs_humanCategoryNullCompute_SBP-spatial.sh

# ------------------------------------------------------
# Compute CFPRs for ensembles (under conventional GO enrichment):
# ------------------------------------------------------
qsub pbs_humanUniformRandomShuffleUniformRandom.sh
qsub pbs_humanUniformRandom.sh
qsub pbs_humanSpatialLag.sh

# Intra-category correlation (but randomized versions are not used):
qsub pbs_humanIntraCorr_raw.sh
