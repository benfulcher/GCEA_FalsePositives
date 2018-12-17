
whatGOIDs = [7215,32612];
whatSpecies = 'mouse';
whatSurrogate = 'spatialLag';
numNullSamples = 2000;
whatCorr = 'Spearman';

categoryScores = struct();
categoryLabels = struct();
[categoryScores.spatialLag,categoryLabels.spatialLag] = CompareNulls(whatGOIDs,whatSpecies,whatSurrogate,numNullSamples);
[categoryScores.randomMap,categoryLabels.randomMap] = CompareNulls(whatGOIDs,whatSpecies,whatSurrogate,numNullSamples);


extraParams = struct();
extraParams.theColors = {[32,178,170]/255,[184,134,11]/255}; % [220,220,220]/255}; % [119,136,153]/255,  [184,134,11]/255, [119,136,153]/255, [220,220,220]/255};
extraParams.customSpot = '';
extraParams.offsetRange = 0.7;
BF_JitteredParallelScatter(categoryScores,true,true,true,extraParams);
ax = gca();
ax.XTick = 1:2;
ax.XTickLabel = categoryLabels;
ylabel('Mean category score')
xlabel('GO category')
title(sprintf('%u nulls of %s',numNullSamples,whatSurrogate))
