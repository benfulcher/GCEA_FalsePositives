function CaseStudyResults(whatSpecies,structFilter)
% Case study of degree enrichment in mouse brain, mouse cortex, and human cortex
%-------------------------------------------------------------------------------

if nargin < 1
    % whatSpecies = 'mouse';
    whatSpecies = 'human';
end
if nargin < 2
    structFilter = 'cortex'; % 'cortex', 'all'
end
corrType = 'Spearman';
%-------------------------------------------------------------------------------

% Set general parameters common to all analyses:
params = GiveMeDefaultParams(whatSpecies,structFilter);

%-------------------------------------------------------------------------------
% Compute results:
%-------------------------------------------------------------------------------
resultTablesDegree = struct();
% (i) random-gene null:
resultTablesDegree.randomGeneNull = NodeSimpleEnrichment(params,'degree',corrType);
% (ii) random phenotype null:
resultTablesDegree.randomMap = PerformEnrichment(params,'degree','randomMap');
% (iii) spatial lag null:
resultTablesDegree.spatialLag = PerformEnrichment(params,'degree','spatialLag');

% List significant categories under each null:
whatPField = 'pValZCorr';
countMe = @(x)sum(resultTablesDegree.(x).(whatPField) < params.e.sigThresh);
fprintf(1,'%u categories significant (%s) for random gene null\n',countMe('randomGeneNull'),whatPField);
fprintf(1,'%u categories significant (%s) for random phenotype null\n',countMe('randomMap'),whatPField);
fprintf(1,'%u categories significant (%s) for spatial-lag null\n',countMe('spatialLag'),whatPField);

%-------------------------------------------------------------------------------
% Assemble a joint table:
%-------------------------------------------------------------------------------
[commonGOIDs,ia,ib] = intersect(resultTablesDegree.randomGeneNull.GOID,resultTablesDegree.randomMap.GOID);

GOName = resultTablesDegree.randomGeneNull.GOName(ia);
GOIDlabel = resultTablesDegree.randomGeneNull.GOIDlabel(ia);
GOID = commonGOIDs;
pValZCorrRandomGene = resultTablesDegree.randomGeneNull.pValZCorr(ia);
pValZCorrRandomMap = resultTablesDegree.randomMap.pValZCorr(ib);
pValZCorrSpatialLag = resultTablesDegree.spatialLag.pValZCorr(ib);

newTable = table(GOName,GOIDlabel,GOID,...
                pValZCorrRandomGene,pValZCorrRandomMap,pValZCorrSpatialLag);
meanScoreSum = newTable.pValZCorrRandomGene + newTable.pValZCorrRandomMap;
[~,ix] = sort(meanScoreSum,'ascend');
newTable = newTable(ix,:);

%-------------------------------------------------------------------------------
% Extra analysis-specific analyses
%-------------------------------------------------------------------------------
switch structFilter
case 'all'
    %-------------------------------------------------------------------------------
    % Do scores correlate with intracategory coexpression?
    %-------------------------------------------------------------------------------
    % Obtain the GOTable for ranksum expression differences in isocortex:
    GOTable_isocortex = NodeSimpleEnrichment(params,'isocortex','all');
    PlotGOScoreScatter(resultTablesDegree.randomGeneNull,GOTable_isocortex,{'meanScore','meanScore'});
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
        fprintf(1,'\n%u/%u: %s\n',i,topWhat,resultTablesDegree.randomGeneNull.GOName{i});
        fprintf(1,'Random-gene: %s = %.2g\n',whatPField,resultTablesDegree.randomGeneNull.(whatPField)(i));
        isHere = find(resultTablesDegree.randomMap.GOID==resultTablesDegree.randomGeneNull.GOID(i));
        fprintf(1,'Random phenotype: %s = %.2g\n',whatPField,resultTablesDegree.randomMap.(whatPField)(isHere));
        isHere = find(resultTablesDegree.randomMap.GOID==resultTablesDegree.spatialLag.GOID(i));
        fprintf(1,'spatialLag: %s = %.2g\n',whatPField,resultTablesDegree.spatialLag.(whatPField)(isHere));
    end
end


f = figure('color','w');
histogram(gScores)
fprintf(1,'Correlations range from %.2f--%.2f\n',min(gScores),max(gScores));

% (ii) spatial null:
shuffleWhat = 'all'; % will be just cortex by definition; given inclusion criterion
resultTablesDegree.mouse_ctx_spatialNull = NodeShuffleEnrichment('degree',...
                        shuffleWhat,numNulls,params.c.structFilter,params);

%===============================================================================


%-------------------------------------------------------------------------------
% How correlated are whole-brain category scores between random gene and spatial nulls?
%-------------------------------------------------------------------------------
PlotGOScoreScatter(resultTablesDegree.mouse_all_randomGeneNull,...
                    resultTablesDegree.mouse_all_spatialNullAll,{'pValCorr','pValZCorr'});
xlabel('degree-corr-randomGeneNull')
ylabel('degree-corr-spatial-null')


%-------------------------------------------------------------------------------
% Results when permuting separately within and between cortical areas
%-------------------------------------------------------------------------------
pVals = resultTablesDegree.mouse_all_spatialNull_twoIsocortex.pValZCorr;
fprintf(1,'p-vals for two-part isocortex constrained null are all > %.3f [%u-nulls]\n',...
                min(pVals),numNulls);
