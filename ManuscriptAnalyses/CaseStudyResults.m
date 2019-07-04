function CaseStudyResults(whatSpecies,structFilter)
% Case study of degree enrichment in mouse brain, mouse cortex, and human cortex
%-------------------------------------------------------------------------------

if nargin < 1
    % whatSpecies = 'mouse';
    whatSpecies = 'human';
end
if nargin < 2
    structFilter = 'all'; % 'cortex', 'all'
end
corrType = 'Spearman';
%-------------------------------------------------------------------------------

% Set general parameters common to all analyses:
params = GiveMeDefaultParams(whatSpecies,structFilter);

%-------------------------------------------------------------------------------
% Compute results:
%-------------------------------------------------------------------------------
resultsTablesDegree = struct();
% (i) random-gene null:
resultsTablesDegree.randomGeneNull = NodeSimpleEnrichment(params,'degree',corrType);
% (ii) random phenotype null:
resultsTablesDegree.randomMap = PerformEnrichment(params,'degree','randomMap');
% (iii) spatial lag null:
resultsTablesDegree.spatialLag = PerformEnrichment(params,'degree','spatialLag');

% List significant categories under each null:
whatPField = 'pValZCorr';
countMe = @(x)sum(resultsTablesDegree.(x).(whatPField) < params.e.sigThresh);
fprintf(1,'%u categories significant (%s) for random gene null\n',countMe('randomGeneNull'),whatPField);
fprintf(1,'%u categories significant (%s) for random phenotype null\n',countMe('randomMap'),whatPField);
fprintf(1,'%u categories significant (%s) for spatial-lag null\n',countMe('spatialLag'),whatPField);

%-------------------------------------------------------------------------------
% Extra analysis-specific analyses
%-------------------------------------------------------------------------------
switch whatAnalysis
case 'wholeBrain'
    %-------------------------------------------------------------------------------
    % Do scores correlate with intracategory coexpression?
    %-------------------------------------------------------------------------------
    % Obtain the GOTable for ranksum expression differences in isocortex:
    GOTable_isocortex = NodeSimpleEnrichment(params,'isocortex','all');
    PlotGOScoreScatter(resultsTablesDegree.randomGeneNull,GOTable_isocortex,{'meanScore','meanScore'});
    xlabel('GO category score (mean Spearman correlation with degree)')
    ylabel('GO category score (-log10 p-value ranksum test isocortex)')

    % Does degree differ cortical/non-cortical?
    [k,structInfo] = ComputeDegree(params.humanOrMouse,true);
    kCortex = k(strcmp(structInfo.divisionLabel,'Isocortex'));
    kNotCotex = k(~strcmp(structInfo.divisionLabel,'Isocortex'));
    fprintf(1,'Cortex: <k> = %.2f, s_k = %.2f\n',mean(kCortex),std(kCortex));
    fprintf(1,'Not-cortex: <k> = %.2f, s_k = %.2f\n',mean(kNotCotex),std(kNotCotex));
    p = ranksum(kCortex,kNotCotex);

case 'cortex'
    % Rearrange to look at p-values for each of the top categories
    topWhat = 20;
    for i = 1:topWhat
        fprintf(1,'\n%u/%u: %s\n',i,topWhat,resultsTablesDegree.randomGeneNull.GOName{i});
        fprintf(1,'Random-gene: %s = %.2g\n',whatPField,resultsTablesDegree.randomGeneNull.(whatPField)(i));
        isHere = find(resultsTablesDegree.randomMap.GOID==resultsTablesDegree.randomGeneNull.GOID(i));
        fprintf(1,'Random phenotype: %s = %.2g\n',whatPField,resultsTablesDegree.randomMap.(whatPField)(isHere));
        isHere = find(resultsTablesDegree.randomMap.GOID==resultsTablesDegree.spatialLag.GOID(i));
        fprintf(1,'spatialLag: %s = %.2g\n',whatPField,resultsTablesDegree.spatialLag.(whatPField)(isHere));
    end
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
