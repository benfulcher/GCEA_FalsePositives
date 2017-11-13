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
whatSpecies = 'mouse';
params = GiveMeDefaultParams(whatSpecies);
params.g.normalizationGene = 'none';
% MEAN:
[resultsTablesMean.mouse_all] = NodeSimpleEnrichment('meanExpression','all','',whatSpecies,params);
[resultsTablesMean.mouse_ctx] = NodeSimpleEnrichment('meanExpression','cortex','',whatSpecies,params);
% VARIANCE:
[resultsTablesVar.mouse_all] = NodeSimpleEnrichment('varExpression','all','',whatSpecies,params);
[resultsTablesVar.mouse_ctx] = NodeSimpleEnrichment('varExpression','cortex','',whatSpecies,params);

%-------------------------------------------------------------------------------
% Then we do human stuff:
whatSpecies = 'human';
params = GiveMeDefaultParams(whatSpecies);
params.g.whatParcellation = 'HCP';
params.g.normalizationInternal = 'none';
params.g.normalizationGene = 'none';
params.g.normalizationRegion = 'none';
% MEAN:
[resultsTablesMean.human_HCP] = NodeSimpleEnrichment('meanExpression','cortex','',whatSpecies,params);
% VARIANCE:
[resultsTablesVar.human_HCP] = NodeSimpleEnrichment('varExpression','cortex','',whatSpecies,params);

%-------------------------------------------------------------------------------
% Visualize results
thresholdSig = 0.05;
% MEAN:
PlotEnrichmentTables(resultsTablesMean,thresholdSig);
title('Enrichment by mean expression across the brain')
% VAR:
PlotEnrichmentTables(resultsTablesVar,thresholdSig);
title('Enrichment by expression variance across the brain')

%===============================================================================
% Cortical differences:
%===============================================================================
resultsTablesCortexDiff = struct();
params = GiveMeDefaultParams('mouse');
[resultsTablesCortexDiff.mouse] = NodeSimpleEnrichment('isocortex','all','','mouse',params);
params = GiveMeDefaultParams('human');
params.g.whatParcellation = 'APARC';
params.g.normalizationInternal = 'none';
params.g.normalizationRegion = 'none';
params.g.normalizationGene = 'zscore';
[resultsTablesCortexDiff.human_APARC,gScore] = NodeSimpleEnrichment('isocortex','all','','human',params);

thresholdSig = 0.5;
PlotEnrichmentTables(resultsTablesCortexDiff,thresholdSig);

%===============================================================================
% PCs of expression variation:
%===============================================================================
resultsTablesPC1 = struct();

% (Note that z-score normalization of columns occurs subsequently within the
% enrichment code)

whatSpecies = 'mouse';
params = GiveMeDefaultParams(whatSpecies);
params.g.normalizationGene = 'mixedSigmoid'; % 'none', 'mixedSigmoid'
params.g.normalizationRegion = 'none'; % 'none', 'mixedSigmoid'
structFilter = 'all';
[resultsTablesPC1.mouse_all] = NodeSimpleEnrichment('genePC',structFilter,'',whatSpecies,params);
structFilter = 'isocortex';
[resultsTablesPC1.mouse_ctx] = NodeSimpleEnrichment('genePC',structFilter,'',whatSpecies,params);

whatSpecies = 'human';
params = GiveMeDefaultParams(whatSpecies);
params.g.whatParcellation = 'HCP';
params.g.normalizationInternal = 'robustSigmoid';
params.g.normalizationGene = 'none';
params.g.normalizationRegion = 'none';
[resultsTablesPC1.human_HCP,gScore] = NodeSimpleEnrichment('genePC','all','',whatSpecies,params);

% Give summary to screen:
thresholdSig = 0.05;
countMe = @(x)sum(resultsTablesPC1.(x).pValCorr < thresholdSig);
fprintf(1,'Found %u (mouse-all), %u (mouse-ctx), %u (human-HCP)\n',...
                countMe('mouse_all'),...
                countMe('mouse_ctx'),...
                countMe('human_HCP'));

% PLOT:
thresholdSig = 0.1;
PlotEnrichmentTables(resultsTablesPC1,thresholdSig);
