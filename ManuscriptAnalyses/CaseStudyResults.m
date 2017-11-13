% The case study on degree-enrichment
%-------------------------------------------------------------------------------

resultsTablesDegree = struct();

% First is degree correlations across the mouse brain:

whatSpecies = 'mouse';
params = GiveMeDefaultParams(whatSpecies);
params.g.normalizationGene = 'none';
corrType = 'Spearman';

% Across the whole brain:
structFilter = 'all';
resultsTablesDegree.mouse_all = NodeSimpleEnrichment('degree',structFilter,corrType,whatSpecies,params);

% Across the cortex only:
structFilter = 'isocortex';
resultsTablesDegree.mouse_ctx = NodeSimpleEnrichment('degree',structFilter,corrType,whatSpecies,params);

thresholdSig = 0.05;
countMe = @(x)sum(resultsTablesDegree.(x).pValCorr < thresholdSig);
fprintf(1,'%u categories significant for whole brain, %u for isocortex\n',...
                        countMe('mouse_all'),countMe('mouse_ctx'));
