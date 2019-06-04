% Case studies of degree enrichment
%-------------------------------------------------------------------------------

resultsTablesDegree = struct();

% Set general parameters common to all analyses:
whatSpecies = 'mouse';
params = GiveMeDefaultParams(whatSpecies);
corrType = 'Spearman';

% Num nulls for shuffle-based enrichment:
numNulls = 20000;
% Category significance threshold for listing out to screen:
thresholdSig = 0.05;

%===============================================================================
% ---1---Across the whole brain:
%===============================================================================
params.c.structFilter = 'all';
% (i) random gene null:
resultsTablesDegree.mouse_randomGeneNull = NodeSimpleEnrichment('degree',...
                                        params.c.structFilter,corrType,params);
% (ii) random phenotype null:
resultsTablesDegree.mouse_randomMap = PerformEnrichment('degree','mouse','randomMap');
% (iii) spatial lag null:
resultsTablesDegree.mouse_spatialLag = PerformEnrichment('degree','mouse','spatialLag');

% Statistics on significance
countMe = @(x)sum(resultsTablesDegree.(x).pValZCorr < thresholdSig);
fprintf(1,'%u categories significant for random gene null\n',countMe('mouse_randomGeneNull'));
fprintf(1,'%u categories significant for random phenotype null\n',countMe('mouse_randomMap'));
fprintf(1,'%u categories significant for spatial-lag null\n',countMe('mouse_spatialLag'));

% Do scores correlate with intracategory coexpression?
% resultsTablesDegree.mouse_all_randomGeneNull.meanScore



% shuffleWhat = 'all'; % random shuffling
% resultsTablesDegree.mouse_all_spatialNullAll = NodeShuffleEnrichment('degree',...
%                 shuffleWhat,numNulls,params.c.structFilter,params);

% (iii) spatial null (cortex constrained):
% shuffleWhat = 'twoIsocortex';
% resultsTablesDegree.mouse_all_spatialNull_twoIsocortex = NodeShuffleEnrichment('degree',...
%                 shuffleWhat,numNulls,params.c.structFilter,params);

%===============================================================================

%-------------------------------------------------------------------------------
% How correlated are degree scores with cortical scores
%-------------------------------------------------------------------------------
resultsTablesCortex = NodeSimpleEnrichment('isocortex',params.c.structFilter,...
                            corrType,params);

% With random gene null scores:
PlotGOScoreScatter(resultsTablesCortex,resultsTablesDegree.mouse_randomGeneNull,{'meanScore','meanScore'});
xlabel('cortex-noncortex')
ylabel('degree-corr-randomgeneNull')


%===============================================================================
% ---2---Across the cortex only:
%===============================================================================
params.c.structFilter = 'isocortex';

% (i) randomGene null:
[resultsTablesDegree.mouse_ctx_randomGene,gScores] = NodeSimpleEnrichment('degree',...
                                        params.c.structFilter,corrType,params);


f = figure('color','w');
histogram(gScores)
fprintf(1,'Correlations range from %.2f--%.2f\n',min(gScores),max(gScores));

% (ii) spatial null:
shuffleWhat = 'all'; % will be just cortex by definition; given inclusion criterion
resultsTablesDegree.mouse_ctx_spatialNull = NodeShuffleEnrichment('degree',...
                        shuffleWhat,numNulls,params.c.structFilter,params);

%===============================================================================


%-------------------------------------------------------------------------------
% How correlated are whole-brain category scores between random gene and spatial nulls?
%-------------------------------------------------------------------------------
PlotGOScoreScatter(resultsTablesDegree.mouse_all_randomGeneNull,...
                    resultsTablesDegree.mouse_all_spatialNullAll,{'pValCorr','pValZCorr'});
xlabel('degree-corr-randomGeneNull')
ylabel('degree-corr-spatial-null')



%-------------------------------------------------------------------------------
% Results when permuting separately within and between cortical areas
%-------------------------------------------------------------------------------
pVals = resultsTablesDegree.mouse_all_spatialNull_twoIsocortex.pValZCorr;
fprintf(1,'p-vals for two-part isocortex constrained null are all > %.3f [%u-nulls]\n',...
                min(pVals),numNulls);
