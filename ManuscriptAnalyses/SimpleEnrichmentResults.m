function SimpleEnrichmentResults(whatAnalysis)
%===============================================================================
% Confound characterization analyses
%===============================================================================

if nargin < 1
    whatAnalysis = 'mean';
end

results = struct();
switch whatAnalysis
case 'mean'
    %===============================================================================
    % Mean expression levels:
    %===============================================================================
    params = GiveMeDefaultParams('mouse');
    results.mouseBrain = NodeSimpleEnrichment(params,'meanExpression');
    params = GiveMeDefaultParams('mouse','cortex');
    results.mouseCtx = NodeSimpleEnrichment(params,'meanExpression');
    params = GiveMeDefaultParams('human');
    results.human = NodeSimpleEnrichment(params,'meanExpression');
case 'var'
    params = GiveMeDefaultParams('mouse','all');
    results.mouseBrain = NodeSimpleEnrichment(params,'varExpression');
    params = GiveMeDefaultParams('mouse','cortex');
    results.mouseCtx = NodeSimpleEnrichment(params,'varExpression');
    params = GiveMeDefaultParams('human');
    results.human = NodeSimpleEnrichment(params,'varExpression');
case 'gradient'
    % (Note that z-score normalization of columns occurs subsequently within the
    % enrichment code)
    params = GiveMeDefaultParams('mouse','all');
    params.g.normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
    params.g.normalizationRegion = 'none'; % 'none', 'mixedSigmoid'
    results.mouseBrain = NodeSimpleEnrichment(params,'genePC');
    params = GiveMeDefaultParams('mouse','cortex');
    params.g.normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
    params.g.normalizationRegion = 'none'; % 'none', 'mixedSigmoid'
    results.mouseCtx = NodeSimpleEnrichment(params,'genePC');
    params = GiveMeDefaultParams('human');
    params.g.normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
    params.g.normalizationRegion = 'none'; % 'none', 'mixedSigmoid'
    results.human = NodeSimpleEnrichment(params,'genePC');
end

[commonGOIDs,ia,ib] = intersect(results.mouseBrain.GOID,results.human.GOID);

GOName = results.mouseBrain.GOName(ia);
GOIDlabel = results.mouseBrain.GOIDlabel(ia);
GOID = commonGOIDs;
sizeMouse = results.mouseBrain.size(ia);
sizeHuman = results.human.size(ib);
meanScoreMouseBrain = results.mouseBrain.meanScore(ia);
meanScoreMouseCortex = results.mouseCtx.meanScore(ia);
meanScoreHuman = results.human.meanScore(ib);
pValZCorrMouseBrain = results.mouseBrain.pValZCorr(ia);
pValZCorrMouseCortex = results.mouseCtx.pValZCorr(ia);
pValZCorrHuman = results.human.pValZCorr(ib);

newTable = table(GOName,GOIDlabel,GOID,sizeMouse,sizeHuman,...
                meanScoreMouseBrain,meanScoreMouseCortex,meanScoreHuman,...
                pValZCorrMouseBrain,pValZCorrMouseCortex,pValZCorrHuman);
meanScoreSum = newTable.pValZCorrMouseBrain + newTable.pValZCorrHuman;
[~,ix] = sort(meanScoreSum,'ascend');
newTable = newTable(ix,:);



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
% PCs of expression variation:
%===============================================================================
resultsTablesPC1 = struct();



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
