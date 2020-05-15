function DistanceConfoundResults()
%-------------------------------------------------------------------------------
% Investigate the enrichment signatures of distance-related confounds
%-------------------------------------------------------------------------------

% ^*^*^*^ Make sure results are computed first using ComputeSpatialEmbeddingScores

% Store results tables in this struct:
results = struct();

%-------------------------------------------------------------------------------
% NB: To make self-correlation make sense, expression needs to be normalized
%-------------------------------------------------------------------------------
% Mouse brain
load(GiveMeDistanceScoreFileName(GiveMeDefaultParams('mouse','all')),'GOTable');
results.mouseBrain = GOTable;
% Mouse cortex:
load(GiveMeDistanceScoreFileName(GiveMeDefaultParams('mouse','cortex')),'GOTable');
results.mouseCtx = GOTable;
% Human
load(GiveMeDistanceScoreFileName(GiveMeDefaultParams('human','cortex')),'GOTable');
results.human = GOTable;
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Assemble a joint table:
[commonGOIDs,ia,ib] = intersect(results.mouseBrain.GOID,results.human.GOID);

GOName = results.mouseBrain.GOName(ia);
GOIDlabel = results.mouseBrain.GOIDlabel(ia);
GOID = commonGOIDs;
meanScoreMouseBrain = results.mouseBrain.meanScore(ia);
pValZCorrMouseBrain = results.mouseBrain.pValZCorr(ia);
meanScoreMouseCtx = results.mouseCtx.meanScore(ia);
pValZCorrMouseCtx = results.mouseCtx.pValZCorr(ia);
meanScoreHuman = results.human.meanScore(ib);
pValZCorrHuman = results.human.pValZCorr(ib);

newTable = table(GOName,GOIDlabel,GOID,...
                meanScoreMouseBrain,meanScoreMouseCtx,meanScoreHuman,...
                pValZCorrMouseBrain,pValZCorrMouseCtx,pValZCorrHuman);
meanScoreSum = newTable.meanScoreMouseBrain + newTable.meanScoreMouseCtx + newTable.meanScoreHuman;
[~,ix] = sort(meanScoreSum,'descend');
newTable = newTable(ix,:);

%-------------------------------------------------------------------------------
% Some basic stats:
fprintf(1,'%.2f%% have a positive autocorrelation score (mouse brain)\n',...
                    mean(meanScoreMouseBrain > 0)*100);
fprintf(1,'%.2f%% have a positive autocorrelation score (human cortex)\n',...
                    mean(meanScoreHuman > 0)*100);

%-------------------------------------------------------------------------------
% Save it to csv file:
IDLabel = newTable.GOIDlabel;
CategoryName = newTable.GOName;
ID = newTable.GOID;
MeanSpatialAutocorrelationScore_MouseBrain = newTable.meanScoreMouseBrain;
MeanSpatialAutocorrelationScore_MouseCortex = newTable.meanScoreMouseCtx;
MeanSpatialAutocorrelationScore_Human = newTable.meanScoreHuman;
% FDRpVal_MouseBrain = newTable.pValZCorrMouseBrain;
% FDRpVal_MouseCortex = newTable.pValZCorrMouseCtx;
% FDRpVal_Human = newTable.pValZCorrHuman;
T = table(CategoryName,IDLabel,ID,MeanSpatialAutocorrelationScore_MouseBrain,...
                        MeanSpatialAutocorrelationScore_MouseCortex,...
                        MeanSpatialAutocorrelationScore_Human);
fileOut = fullfile('SupplementaryTables','SpatialAutocorrelationScores.csv');
writetable(T,fileOut,'Delimiter',',','QuoteStrings',true);
fprintf(1,'Saved spatial autocorrelation scores to %s\n',fileOut);

%-------------------------------------------------------------------------------
% Visualize any overlapping spatial signatures:
% thresholdSig = [0.05,2,2];
% PlotEnrichmentTables(results,thresholdSig);

%-------------------------------------------------------------------------------
% Visualize specific categories:
% params = GiveMeDefaultParams('mouse');
% whatCategoryIndex = 1; % (NB: index not ID)
% VisualizeDistanceEnrichment(results.mouseBrain,whatCategoryIndex,params);
