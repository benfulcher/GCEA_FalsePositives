% DistanceConfoundResults Investigate the enrichment signatures of distance-related confounds

% Store results tables in this struct:
results = struct();

%-------------------------------------------------------------------------------
% Mouse
%-------------------------------------------------------------------------------
params = GiveMeDefaultParams('mouse');
% To make GCC scores make sense -- expression needs to be [0,1] normalized:
params.g.normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
params.g.normalizationRegion = 'zscore'; % 'none', 'zscore'

% ---Whole brain:
params.c.structFilter = 'all';
results.mouse_all = geneEnrichmentDistance(params);
% ---Just isocortex:
params.c.structFilter = 'isocortex';
results.mouse_ctx = geneEnrichmentDistance(params);
% ---Just non-cortical areas:
params.c.structFilter = 'notCortex';
results.mouse_notCtx = geneEnrichmentDistance(params);

%-------------------------------------------------------------------------------
% Investigate specific categories through specific visualizations:
whatCategoryIndex = 1; % (NB: index not ID)
VisualizeDistanceEnrichment(results.mouse_all,whatCategoryIndex,params);


%-------------------------------------------------------------------------------
% Human
%-------------------------------------------------------------------------------
params = GiveMeDefaultParams('human');
params.g.normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
params.g.normalizationRegion = 'zscore'; % 'none', 'zscore'

% ---Whole brain:
params.c.structFilter = 'all';
results.human_all = geneEnrichmentDistance(params);

% ---Just cortex:
params.c.structFilter = 'cortex';
results.human_ctx = geneEnrichmentDistance(params);

%-------------------------------------------------------------------------------
% Visualize together:
thresholdSig = [0.05,2,3];
PlotEnrichmentTables(results,thresholdSig);
title('Enrichment by expression variance across the brain')
cB = colorbar();
