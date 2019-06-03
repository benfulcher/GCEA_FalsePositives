% PlotCategoryNullCompare
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Parameters:
whatGOIDs = [7215,32612];
whatSpecies = 'mouse';
numNullSamples = 2000;
whatCorr = 'Spearman';

%-------------------------------------------------------------------------------
% Plot information about a category
categoryWhat = 2;
params = GiveMeDefaultParams(whatSpecies);
PlotCategoryIntraCorr(whatGOIDs(categoryWhat),params,whatCorr);
PlotCategoryExpression(whatGOIDs(categoryWhat),params);

%-------------------------------------------------------------------------------

categoryScores = struct();
categoryLabels = struct();
fprintf(1,'Computing distribution of null spatial-lag category scores for %u GO categories\n',length(whatGOIDs));
[categoryScores.spatialLag,categoryLabels.spatialLag] = CompareNulls(whatGOIDs,whatSpecies,'spatialLag',numNullSamples,true);
fprintf(1,'Computing distribution of null random-map category scores for %u GO categories\n',length(whatGOIDs));
[categoryScores.randomMap,categoryLabels.randomMap] = CompareNulls(whatGOIDs,whatSpecies,'randomMap',numNullSamples,false);

%-------------------------------------------------------------------------------
categoryScoresTogether = [categoryScores.spatialLag; categoryScores.randomMap];
categoryLabelsTogether = [categoryLabels.spatialLag; categoryLabels.randomMap];

% Tests of variance:
[h,p] = vartest2(categoryScores.spatialLag{1},categoryScores.spatialLag{2})
[h,p] = vartest2(categoryScores.randomMap{1},categoryScores.randomMap{2})

%-------------------------------------------------------------------------------
% Violin plots:
f = figure('color','w');
extraParams = struct();
extraParams.theColors = {[32,178,170]/255,[184,134,11]/255,brighten([32,178,170]/255,+0.5),brighten([184,134,11]/255,+0.5)}; % [220,220,220]/255}; % [119,136,153]/255, [119,136,153]/255, [220,220,220]/255};
extraParams.customSpot = '';
extraParams.offsetRange = 0.7;
BF_JitteredParallelScatter(categoryScoresTogether,true,true,false,extraParams);
ax = gca();
ax.XTick = 1:4;
ax.XTickLabel = categoryLabelsTogether;
ylabel('Mean category score')
xlabel('GO category')
title(sprintf('%u nulls',numNullSamples))
