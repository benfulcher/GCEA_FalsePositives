%===============================================================================
% Confound characterization analyses
%===============================================================================

%===============================================================================
% Mean/variance expression levels:
%===============================================================================
resultsTablesMean = struct();
resultsTablesVar = struct();

%-------------------------------------------------------------------------------
% First we do mouse stuff:
params = GiveMeDefaultParams('mouse');
params.g.normalizationGene = 'none';
% MEAN:
[resultsTablesMean.mouse_all] = NodeSimpleEnrichment('meanExpression','all','','mouse',params);
[resultsTablesMean.mouse_ctx] = NodeSimpleEnrichment('meanExpression','cortex','','mouse',params);
% VARIANCE:
[resultsTablesVar.mouse_all] = NodeSimpleEnrichment('varExpression','all','','mouse',params);
[resultsTablesVar.mouse_ctx] = NodeSimpleEnrichment('varExpression','cortex','','mouse',params);

%-------------------------------------------------------------------------------
% Then we do human stuff:
params = GiveMeDefaultParams('human');
params.g.whatParcellation = 'HCP';
params.g.normalizationInternal = 'none';
params.g.normalizationGene = 'none';
params.g.normalizationRegion = 'none';
% MEAN:
[resultsTablesMean.human_HCP] = NodeSimpleEnrichment('meanExpression','cortex','','human',params);
% VARIANCE:
[resultsTablesVar.human_HCP] = NodeSimpleEnrichment('varExpression','cortex','','human',params);

%-------------------------------------------------------------------------------
% Visualize results
% MEAN:
PlotEnrichmentTables(resultsTablesMean,0.05);
title('Enrichment by mean expression across the brain')
% VAR:
PlotEnrichmentTables(resultsTablesVar,0.05);
title('Enrichment by expression variance across the brain')

%===============================================================================
% Spatial variation
%===============================================================================
resultsTablesCortexDiff = struct();
params = GiveMeDefaultParams('mouse');
[resultsTablesCortexDiff.mouse] = NodeSimpleEnrichment('isocortex','all','','mouse',params);
params = GiveMeDefaultParams('human');
params.g.whatParcellation = 'APARC';
params.g.normalizationInternal = 'robustSigmoid';
params.g.normalizationRegion = 'none';
params.g.normalizationGene = 'zscore';
[resultsTablesCortexDiff.human_APARC] = NodeSimpleEnrichment('isocortex','all','','human',params);
