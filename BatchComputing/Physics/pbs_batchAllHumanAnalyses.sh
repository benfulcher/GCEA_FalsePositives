#!/bin/bash

# ------------------------------------------------------
# Compute null distributions under different phenotype ensembles:
# ------------------------------------------------------
qsub pbs_humanCategoryNullCompute_SBP-rand.sh
qsub pbs_humanCategoryNullCompute_SBP-spatial.sh

# ------------------------------------------------------
# Compute CFPRs for ensembles (under conventional GO enrichment):
# ------------------------------------------------------
qsub pbs_humanCFPR-Ref.sh
qsub pbs_humanCFPR_SBP-rand.sh
qsub pbs_humanSpatialLag.sh

# Intra-category correlation (but randomized versions are not used):
qsub pbs_humanIntraCorr_raw.sh
