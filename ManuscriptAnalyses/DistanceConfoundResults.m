%-------------------------------------------------------------------------------
% DistanceConfoundResults
%-------------------------------------------------------------------------------
% Investigate the enrichment signatures of distance-related confounds
%-------------------------------------------------------------------------------

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
results.mouseBrain = geneEnrichmentDistance(params);

% ---Just isocortex:
params.c.structFilter = 'isocortex';
results.mouseCtx = geneEnrichmentDistance(params);

% ---Just non-cortical areas:
% params.c.structFilter = 'notCortex';
% results.mouse_notCtx = geneEnrichmentDistance(params);

%-------------------------------------------------------------------------------
% Human
%-------------------------------------------------------------------------------
params = GiveMeDefaultParams('human');
params.g.normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
params.g.normalizationRegion = 'zscore'; % 'none', 'zscore'

% ---Just cortex:
results.human = geneEnrichmentDistance(params);


%-------------------------------------------------------------------------------
% Visualize any overlapping spatial signatures:
thresholdSig = [0.05,2,2];
PlotEnrichmentTables(results,thresholdSig);

%-------------------------------------------------------------------------------
% Investigate specific categories through specific visualizations:
whatCategoryIndex = 1; % (NB: index not ID)
VisualizeDistanceEnrichment(results.mouse_all,whatCategoryIndex,params);
