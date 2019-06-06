function CaseStudyResults(whatSpecies,whatAnalysis)
%-------------------------------------------------------------------------------
% Case studies of degree enrichment in mouse brain, mouse cortex, and human cortex
%-------------------------------------------------------------------------------

if nargin < 1
    whatSpecies = 'mouse';
end
if nargin < 2
    whatAnalyais = 'wholeBrain';
end


%-------------------------------------------------------------------------------

% Set general parameters common to all analyses:
corrType = 'Spearman';
thresholdSig = 0.05; % Category significance threshold for listing out to screen:
params = GiveMeDefaultParams(whatSpecies);

switch whatAnalysis
case 'wholeBrain'
    params.c.structFilter = 'all';
    params.g.structFilter = 'all';
    %===============================================================================
    % ---1---Across the whole brain:
    %===============================================================================
    resultsTablesDegree = struct();
    % (i) random-gene null:
    resultsTablesDegree.mouse_randomGeneNull = NodeSimpleEnrichment(params,'degree',corrType);
    % (ii) random phenotype null:
    resultsTablesDegree.mouse_randomMap = PerformEnrichment(params,'degree','randomMap');
    % (iii) spatial lag null:
    resultsTablesDegree.mouse_spatialLag = PerformEnrichment(params,'degree','spatialLag');

    % Statistics on significance
    countMe = @(x)sum(resultsTablesDegree.(x).pValZCorr < thresholdSig);
    fprintf(1,'%u categories significant for random gene null\n',countMe('mouse_randomGeneNull'));
    fprintf(1,'%u categories significant for random phenotype null\n',countMe('mouse_randomMap'));
    fprintf(1,'%u categories significant for spatial-lag null\n',countMe('mouse_spatialLag'));

    %-------------------------------------------------------------------------------
    % Do scores correlate with intracategory coexpression?
    %-------------------------------------------------------------------------------
    % Obtain the GOTable for ranksum expression differences in isocortex:
    GOTable_isocortex = NodeSimpleEnrichment(params,'isocortex','all');
    PlotGOScoreScatter(resultsTablesDegree.mouse_randomGeneNull,GOTable_isocortex,{'meanScore','meanScore'});
    xlabel('GO category score (mean Spearman correlation with degree)')
    ylabel('GO category score (-log10 p-value ranksum test isocortex)')

    % Does degree differ cortical/non-cortical?
    [k,structInfo] = ComputeDegree(params.humanOrMouse,true);
    kCortex = k(strcmp(structInfo.divisionLabel,'Isocortex'));
    kNotCotex = k(~strcmp(structInfo.divisionLabel,'Isocortex'));
    fprintf(1,'Cortex: <k> = %.2f, s_k = %.2f\n',mean(kCortex),std(kCortex));
    fprintf(1,'Not-cortex: <k> = %.2f, s_k = %.2f\n',mean(kNotCotex),std(kNotCotex));
    p = ranksum(kCortex,kNotCotex)


case 'cortexOnly'
    %===============================================================================
    % ---2---Across the cortex only:
    %===============================================================================
    params.c.structFilter = 'cortex';
    params.g.structFilter = 'cortex';

    resultsTablesDegree = struct();
    % (i) random-gene null:
    resultsTablesDegree.mouse_randomGeneNull = NodeSimpleEnrichment(params,'degree',corrType);
    % (ii) random phenotype null:
    resultsTablesDegree.mouse_randomMap = PerformEnrichment(params,'degree','randomMap');
    % (iii) spatial lag null:
    resultsTablesDegree.mouse_spatialLag = PerformEnrichment(params,'degree','spatialLag');

    % Statistics on significance
    countMe = @(x)sum(resultsTablesDegree.(x).pValZCorr < thresholdSig);
    fprintf(1,'%u categories significant for random gene null\n',countMe('mouse_randomGeneNull'));
    fprintf(1,'%u categories significant for random phenotype null\n',countMe('mouse_randomMap'));
    fprintf(1,'%u categories significant for spatial-lag null\n',countMe('mouse_spatialLag'));

end


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
