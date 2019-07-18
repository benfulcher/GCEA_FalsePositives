%-------------------------------------------------------------------------------
% DistanceConfoundResults
%-------------------------------------------------------------------------------
% Investigate the enrichment signatures of distance-related confounds
%-------------------------------------------------------------------------------

% Store results tables in this struct:
results = struct();

% NB: To make self-correlation make sense, expression needs to be normalized

%-------------------------------------------------------------------------------
% Mouse brain
params = GiveMeDefaultParams('mouse','all');
params.g.normalizationGene = 'zscore';
params.g.normalizationRegion = 'zscore';
results.mouseBrain = geneEnrichmentDistance(params);

% Mouse cortex:
params = GiveMeDefaultParams('mouse','cortex');
params.g.normalizationGene = 'zscore';
params.g.normalizationRegion = 'zscore';
results.mouseCtx = geneEnrichmentDistance(params);

% Human
params = GiveMeDefaultParams('human','cortex');
params.g.normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
params.g.normalizationRegion = 'zscore'; % 'none', 'zscore'
results.human = geneEnrichmentDistance(params);

% Mouse non-cortex:
% params.c.structFilter = 'notCortex';
% results.mouse_notCtx = geneEnrichmentDistance(params);

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
% Visualize any overlapping spatial signatures:
thresholdSig = [0.05,2,2];
PlotEnrichmentTables(results,thresholdSig);

%-------------------------------------------------------------------------------
% Investigate specific categories through specific visualizations:
whatCategoryIndex = 1; % (NB: index not ID)
VisualizeDistanceEnrichment(results.mouseBrain,whatCategoryIndex,params);
