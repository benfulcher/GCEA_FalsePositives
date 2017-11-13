% The case study on degree-enrichment
%-------------------------------------------------------------------------------

resultsTablesDegree = struct();

% Set general parameters common to all analyses:
whatSpecies = 'mouse';
params = GiveMeDefaultParams(whatSpecies);
params.g.normalizationGene = 'none';
corrType = 'Spearman';

% Num nulls for shuffle-based enrichment:
numNulls = 20;
% For listing significant categories out to screen:
thresholdSig = 0.05;

%===============================================================================
% ---1---Across the whole brain:
%===============================================================================
params.c.structFilter = 'all';
% (i) random gene null:
resultsTablesDegree.mouse_all_randomGeneNull = NodeSimpleEnrichment('degree',...
                        params.c.structFilter,corrType,whatSpecies,params);

% (ii) spatial null (all):
shuffleWhat = 'all';
resultsTablesDegree.mouse_all_spatialNullAll = NodeShuffleEnrichment('degree',...
                shuffleWhat,numNulls,params.c.structFilter,whatSpecies,params);

% (iii) spatial null (cortex constrained):
shuffleWhat = 'twoIsocortex';
resultsTablesDegree.mouse_all_spatialNull_twoIsocortex = NodeShuffleEnrichment('degree',...
                shuffleWhat,numNulls,params.c.structFilter,whatSpecies,params);

%===============================================================================
countMe = @(x)sum(resultsTablesDegree.(x).pValCorr < thresholdSig);
fprintf(1,'%u categories significant for whole brain, %u for isocortex\n',...
                        countMe('mouse_all'),countMe('mouse_ctx'));

%-------------------------------------------------------------------------------
% How correlated are whole-brain category scores between random gene and spatial nulls?
%-------------------------------------------------------------------------------
PlotGOScoreScatter(resultsTablesDegree.mouse_all_randomGeneNull,...
                    resultsTablesDegree.mouse_all_spatialNullAll);
xlabel('degree-corr-randomGeneNull')
ylabel('degree-corr-spatial-null')

%-------------------------------------------------------------------------------
% How correlated are degree scores with cortical scores
%-------------------------------------------------------------------------------
resultsTablesCortex = NodeSimpleEnrichment('isocortex',params.c.structFilter,...
                            corrType,whatSpecies,params);

% With random gene null scores:
PlotGOScoreScatter(resultsTablesCortex,resultsTablesDegree.mouse_all_randomGeneNull);
xlabel('cortex-noncortex')
ylabel('degree-corr-randomgeneNull')

% With spatial null scores:
PlotGOScoreScatter(resultsTablesCortex,resultsTablesDegree.mouse_all_spatialNullAll);
xlabel('cortex-noncortex')
ylabel('degree-corr-spatialNull')

%-------------------------------------------------------------------------------
% Results when permuting separately within and between cortical areas
%-------------------------------------------------------------------------------
pVals = resultsTablesDegree.mouse_all_spatialNull_twoIsocortex.pValZCorr;
fprintf(1,'p-vals for two-part isocortex constrained null are all > %.3f [%u-nulls]\n',...
                min(pVals),numNulls);

%===============================================================================
% ---2---Across the cortex only:
%===============================================================================
params.c.structFilter = 'isocortex';

% (i) randomGene null:
[resultsTablesDegree.mouse_ctx_randomGene,gScores] = NodeSimpleEnrichment('degree',...
                        params.c.structFilter,corrType,whatSpecies,params);

f = figure('color','w');
histogram(gScores)
fprintf(1,'Correlations range from %.2f--%.2f\n',min(gScores),max(gScores));

% (ii) spatial null:
shuffleWhat = 'all'; % will be just cortex by definition; given inclusion criterion
resultsTablesDegree.mouse_ctx_spatialNull = NodeShuffleEnrichment('degree',...
            shuffleWhat,numNulls,params.c.structFilter,whatSpecies,params);
