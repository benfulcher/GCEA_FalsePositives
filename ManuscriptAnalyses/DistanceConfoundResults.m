% --DistanceConfoundResults--
% Investigate the enrichment signatures of distance-related confounds:

% Store results here:
resultsTablesDist = struct();

%-------------------------------------------------------------------------------
% Set up general GCC compute parameters:
%-------------------------------------------------------------------------------
GCCparams = struct();
GCCparams.whatCorr = 'Spearman'; % 'Pearson', 'Spearman'
GCCparams.pValOrStat = 'stat'; % 'pval','stat'
GCCparams.thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
GCCparams.absType = 'neg';
GCCparams.onlyConnections = false; % only look where there are structural connections
GCCparams.regressDistance = false; % whether to regress distance

%-------------------------------------------------------------------------------
% MOUSE:
%-------------------------------------------------------------------------------
whatSpecies = 'mouse';
% Other parameters:
params = GiveMeDefaultParams(whatSpecies);
% To make GCC scores make sense -- expression needs to be [0,1] normalized:
params.g.normalizationGene = 'mixedSigmoid'; % 'none', 'mixedSigmoid'
params.g.normalizationRegion = 'none'; % 'none', 'zscore'

%----Get the results (for whole-brain and cortex):
structFilter = 'all';
resultsTablesDist.mouse_all = geneEnrichmentDistance(structFilter,whatSpecies,params,GCCparams);
structFilter = 'isocortex';
resultsTablesDist.mouse_ctx = geneEnrichmentDistance(structFilter,whatSpecies,params,GCCparams);

%-------------------------------------------------------------------------------
% Now Human:
whatSpecies = 'human';
params = GiveMeDefaultParams(whatSpecies);
params.g.whatParcellation = 'HCP';
params.g.normalizationInternal = 'robustSigmoid';
params.g.normalizationGene = 'mixedSigmoid'; % 'none', 'mixedSigmoid'
params.g.normalizationRegion = 'none'; % 'none', 'zscore'

resultsTablesDist.human_HCP = geneEnrichmentDistance(structFilter,whatSpecies,params,GCCparams);

%-------------------------------------------------------------------------------
% Visualize together:
thresholdSig = 0.05;
PlotEnrichmentTables(resultsTablesDist,thresholdSig);
title('Enrichment by expression variance across the brain')
