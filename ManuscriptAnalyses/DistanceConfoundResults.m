% DistanceConfoundResults   Investigate the enrichment signatures of distance-related confounds

% Store results tables in this struct:
resultsTablesDist = struct();

%-------------------------------------------------------------------------------
% Mouse:
params = GiveMeDefaultParams('mouse');
% To make GCC scores make sense -- expression needs to be [0,1] normalized:
params.g.normalizationGene = 'mixedSigmoid'; % 'none', 'mixedSigmoid'
params.g.normalizationRegion = 'none'; % 'none', 'zscore'

%----Get the results (for whole-brain and cortex):
params.c.structFilter = 'all';
resultsTablesDist.mouse_all = geneEnrichmentDistance(params);
params.c.structFilter = 'isocortex';
resultsTablesDist.mouse_ctx = geneEnrichmentDistance(params);

%-------------------------------------------------------------------------------
% Human:
params = GiveMeDefaultParams('human');
params.g.whatParcellation = 'HCP';
params.g.normalizationInternal = 'robustSigmoid';
params.g.normalizationGene = 'mixedSigmoid'; % 'none', 'mixedSigmoid'
params.g.normalizationRegion = 'none'; % 'none', 'zscore'

resultsTablesDist.human_HCP = geneEnrichmentDistance(params);

%-------------------------------------------------------------------------------
% Visualize together:
thresholdSig = 0.05;
PlotEnrichmentTables(resultsTablesDist,thresholdSig);
title('Enrichment by expression variance across the brain')
