
%-------------------------------------------------------------------------------
% Parameters:
whatGOIDs = [7215,32612];
whatSpecies = 'mouse';
numNullSamples = 2000;
whatCorr = 'Spearman';

%-------------------------------------------------------------------------------
params = GiveMeDefaultParams(whatSpecies);
PlotCategoryIntraCorr(whatGOIDs(2),params,whatCorr);
PlotCategoryExpression(whatGOIDs(2),params);

%-------------------------------------------------------------------------------

categoryScores = struct();
categoryLabels = struct();
[categoryScores.spatialLag,categoryLabels.spatialLag] = CompareNulls(whatGOIDs,whatSpecies,'spatialLag',numNullSamples);
[categoryScores.randomMap,categoryLabels.randomMap] = CompareNulls(whatGOIDs,whatSpecies,'spatialShuffle',numNullSamples);

%-------------------------------------------------------------------------------
categoryScoresTogether = [categoryScores.spatialLag;categoryScores.randomMap];
categoryLabelsTogether = [categoryLabels.spatialLag;categoryLabels.randomMap];

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
