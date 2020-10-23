function newTable = CaseStudyResults(whatSpecies,structFilter,enrichWhat)
% Compute case study results: e.g., node degree enrichment
%-------------------------------------------------------------------------------
if nargin < 1
    whatSpecies = 'mouse'; % 'mouse', 'human'
end
if nargin < 2
    structFilter = 'all'; % 'cortex', 'all'
end
if nargin < 3
    enrichWhat = 'degree';
end
%-------------------------------------------------------------------------------

% Set general parameters common to all analyses:
params = GiveMeDefaultParams(whatSpecies,structFilter);

% Set up the gene-expression data of interest:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
geneDataStruct = struct();
geneDataStruct.expressionMatrix = geneData;
geneDataStruct.entrezIDs = geneInfo.entrez_id;

%-------------------------------------------------------------------------------
% Get the phenotype of interest and match:
switch enrichWhat
    case 'degree'
        doBinarize = true;
        phenotypeVector = ComputeDegree(params,doBinarize);
    otherwise
        % A cell density map
        [phenotypeVector,ia] = MatchMeCellDensity(structInfo,enrichWhat);
        structInfo = structInfo(ia,:);
        geneDataStruct.expressionMatrix = geneData(ia,:);
end

%-------------------------------------------------------------------------------
% Compute results:
%-------------------------------------------------------------------------------
resultTables = struct();

% -------
% (i) Load results of conventional null (precomputed using gene-score resampling):
% -------
% resultTables.randomGeneNull = NodeSimpleEnrichment(params,'degree',true);
fprintf(1,'Random gene null\n');
load(GiveMeSimpleEnrichmentOutputFile(params,enrichWhat),'GOTable');
resultTables.randomGeneNull = GOTable;
clear('GOTable')

% -------
% (ii) Random phenotype null:
% -------
fprintf(1,'Random phenotype null\n');
% Load the null-phenotype enrichment results (from ComputeAllCategoryNulls):
params.e.whatEnsemble = 'randomMap';
fileNullEnsembleResults = GiveMeEnsembleEnrichmentOutputFileName(params);
% Compute enrichment relative to the null phenotypes:
resultTables.randomMap = EnsembleEnrichment(geneDataStruct,fileNullEnsembleResults,phenotypeVector);
% [geneData,geneInfo,structInfo] = LoadMeG(params.g);
% ListCategories(geneInfo,GOTablePhenotype,20,'pValZ');

% -------
% (iii) Spatial-lag null:
% -------
fprintf(1,'Spatially autocorrelated phenotype null\n');
% Load the null-phenotype enrichment results (from ComputeAllCategoryNulls):
params.e.whatEnsemble = 'customEnsemble';
fileNullEnsembleResults = GiveMeEnsembleEnrichmentOutputFileName(params);
% Compute enrichment relative to the phenotype-based nulls:
resultTables.spatialLag = EnsembleEnrichment(geneDataStruct,fileNullEnsembleResults,phenotypeVector);

%-------------------------------------------------------------------------------
% List significant categories under each null:
whatPField = 'pValZCorr';
countMe = @(x)sum(resultTables.(x).(whatPField) < params.e.sigThresh);
% Can get lost in all the outputs, so make it clear:
fprintf(1,'\n-----------------------------\n');
fprintf(1,'-----------------------------\n');
fprintf(1,'%s-%s-%s\n',whatSpecies,structFilter,enrichWhat);
fprintf(1,'-----------------------------\n');
fprintf(1,'%u categories significant (%s) for random-gene null\n',countMe('randomGeneNull'),whatPField);
fprintf(1,'%u categories significant (%s) for random-phenotype null\n',countMe('randomMap'),whatPField);
fprintf(1,'%u categories significant (%s) for spatial-lag null\n',countMe('spatialLag'),whatPField);
fprintf(1,'-----------------------------\n');
fprintf(1,'-----------------------------\n\n');

%-------------------------------------------------------------------------------
% Assemble a joint table:
%-------------------------------------------------------------------------------
[commonGOIDs,ia,ib] = intersect(resultTables.randomGeneNull.GOID,resultTables.randomMap.GOID);

GOName = resultTables.randomGeneNull.GOName(ia);
GOIDlabel = resultTables.randomGeneNull.GOIDlabel(ia);
GOID = commonGOIDs;
% Raw p-values:
pValZRandomGene = resultTables.randomGeneNull.pValZ(ia);
pValZRandomMap = resultTables.randomMap.pValZ(ib);
pValZSpatialLag = resultTables.spatialLag.pValZ(ib);
% Corrected p-values:
pValZCorrRandomGene = resultTables.randomGeneNull.pValZCorr(ia);
pValZCorrRandomMap = resultTables.randomMap.pValZCorr(ib);
pValZCorrSpatialLag = resultTables.spatialLag.pValZCorr(ib);

newTable = table(GOName,GOIDlabel,GOID,...
                pValZRandomGene,pValZRandomMap,pValZSpatialLag,...
                pValZCorrRandomGene,pValZCorrRandomMap,pValZCorrSpatialLag);
meanScoreSum = newTable.pValZRandomGene + newTable.pValZRandomMap;
[~,ix] = sort(meanScoreSum,'ascend');
newTable = newTable(ix,:);

%-------------------------------------------------------------------------------
% Save out to .csv for a Supplementary Table for the paper:
IDLabel = newTable.GOIDlabel;
CategoryName = newTable.GOName;
ID = newTable.GOID;
FDRpValue_randomGene = newTable.pValZCorrRandomGene;
FDRpValue_SBPrandom = newTable.pValZCorrRandomMap;
FDRpValue_SBPspatial = newTable.pValZCorrSpatialLag;

T = table(CategoryName,IDLabel,ID,FDRpValue_randomGene,FDRpValue_SBPrandom,...
                    FDRpValue_SBPspatial);
fileOut = fullfile('SupplementaryTables',sprintf('EnrichmentThreeWays_%s_%s_%s.csv',whatSpecies,structFilter,enrichWhat));
writetable(T,fileOut,'Delimiter',',','QuoteStrings',true);
fprintf(1,'Saved all CFPR results to %s\n',fileOut);

%-------------------------------------------------------------------------------
% Extra analysis-specific analyses
%-------------------------------------------------------------------------------
switch structFilter
case 'all'
    %-------------------------------------------------------------------------------
    % Do scores correlate with intracategory coexpression?
    %-------------------------------------------------------------------------------
    % Obtain the GOTable for ranksum expression differences in isocortex:
    GOTable_isocortex = NodeSimpleEnrichment(params,'isocortex',true);
    PlotGOScoreScatter(resultTables.randomGeneNull,GOTable_isocortex,{'meanScore','meanScore'});
    xlabel('GO category score (mean Spearman correlation with degree)')
    ylabel('GO category score (-log10 p-value ranksum test isocortex)')

    % Does degree differ cortical/non-cortical?
    [k,structInfo] = ComputeDegree(params,true);
    kCortex = k(strcmp(structInfo.divisionLabel,'Isocortex'));
    kNotCotex = k(~strcmp(structInfo.divisionLabel,'Isocortex'));
    fprintf(1,'Cortex: <k> = %.2f, s_k = %.2f\n',mean(kCortex),std(kCortex));
    fprintf(1,'Not-cortex: <k> = %.2f, s_k = %.2f\n',mean(kNotCotex),std(kNotCotex));
    p = ranksum(kCortex,kNotCotex);

case 'cortex'
    % Rearrange to look at p-values for each of the top categories
    % topWhat = 20;
    % for i = 1:topWhat
    %     fprintf(1,'\n%u/%u: %s\n',i,topWhat,resultTables.randomGeneNull.GOName{i});
    %     fprintf(1,'Random-gene: %s = %.2g\n',whatPField,resultTables.randomGeneNull.(whatPField)(i));
    %     isHere = find(resultTables.randomMap.GOID==resultTables.randomGeneNull.GOID(i));
    %     fprintf(1,'Random phenotype: %s = %.2g\n',whatPField,resultTables.randomMap.(whatPField)(isHere));
    %     isHere = find(resultTables.randomMap.GOID==resultTables.spatialLag.GOID(i));
    %     fprintf(1,'spatialLag: %s = %.2g\n',whatPField,resultTables.spatialLag.(whatPField)(isHere));
    % end
end


% f = figure('color','w');
% histogram(gScores)
% fprintf(1,'Correlations range from %.2f--%.2f\n',min(gScores),max(gScores));
%
% % (ii) spatial null:
% shuffleWhat = 'all'; % will be just cortex by definition; given inclusion criterion
% resultTables.mouse_ctx_spatialNull = NodeShuffleEnrichment('degree',...
%                         shuffleWhat,numNulls,params.c.structFilter,params);
%
% %===============================================================================
%
%
% %-------------------------------------------------------------------------------
% % How correlated are whole-brain category scores between random gene and spatial nulls?
% %-------------------------------------------------------------------------------
% PlotGOScoreScatter(resultTables.mouse_all_randomGeneNull,...
%                     resultTables.mouse_all_spatialNullAll,{'pValCorr','pValZCorr'});
% xlabel('degree-corr-randomGeneNull')
% ylabel('degree-corr-spatial-null')
%
%
% %-------------------------------------------------------------------------------
% % Results when permuting separately within and between cortical areas
% %-------------------------------------------------------------------------------
% pVals = resultTables.mouse_all_spatialNull_twoIsocortex.pValZCorr;
% fprintf(1,'p-vals for two-part isocortex constrained null are all > %.3f [%u-nulls]\n',...
%                 min(pVals),numNulls);

end
