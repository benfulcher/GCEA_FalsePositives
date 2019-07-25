% PlotCategoryNullCompare
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Find two good candidates:

% Load intra-category correlation scores:
whatSpecies = 'mouse';
whatShuffle = 'geneShuffle'; % 'geneShuffle', 'independentSpatialShuffle'
whatIntraStat = 'raw';
numNullSamples_intraCorr = 20000; % (Intra_*_*_20000.mat)
fileNameIntra = sprintf('Intra_%s_%s_%s_%u.mat',whatSpecies,whatShuffle,whatIntraStat,numNullSamples_intraCorr);
resultsIntra = load(fileNameIntra);
categorySizeRange = [40,42];
isInSizeRange = (resultsIntra.resultsTable.size >= categorySizeRange(1)) & (resultsIntra.resultsTable.size <= categorySizeRange(2))
resultsIntraSized = resultsIntra.resultsTable(isInSizeRange,:)

%-------------------------------------------------------------------------------
% Parameters:
% whatGOIDs = [7215,32612];
whatGOIDs = [61001,31638];
numNullSamples = 20000;
whatCorr = 'Spearman';

% Assign colors:
theColors = GiveMeColors('twoGOCategories');
% theColors = [[32,178,170]/255;[184,134,11]/255];
% [220,220,220]/255}; % [119,136,153]/255, [119,136,153]/255, [220,220,220]/255};


%===============================================================================
% Plot information about a category
categoryWhat = 2;
params = GiveMeDefaultParams(whatSpecies);
PlotCategoryIntraCorr(whatGOIDs(categoryWhat),params,whatCorr);
PlotCategoryExpression(whatGOIDs(categoryWhat),params);

%-------------------------------------------------------------------------------
% FPSR for this category?
numNullSamples_surrogate = 10000;
FPSR_random = SurrogateEnrichmentProcess(whatSpecies,numNullSamples_surrogate,'randomUniform','');
FPSR_thisCategory = FPSR_random.sumUnderSig(FPSR_random.GOID==whatGOIDs(categoryWhat));
FPSR_random(4285,:)
%===============================================================================

%===============================================================================
% Plot FPSR distributions
%-------------------------------------------------------------------------------
categoryScores = struct();
categoryLabels = struct();
fprintf(1,'Computing distribution of null spatial-lag category scores for %u GO categories\n',length(whatGOIDs));
[categoryScores.spatialLag,categoryLabels.spatialLag] = CompareNulls(whatGOIDs,whatSpecies,'spatialLag',numNullSamples,false);
fprintf(1,'Computing distribution of null random-map category scores for %u GO categories\n',length(whatGOIDs));
[categoryScores.randomMap,categoryLabels.randomMap] = CompareNulls(whatGOIDs,whatSpecies,'randomMap',numNullSamples,false);

%-------------------------------------------------------------------------------
% Tests of variance:
[h,p] = vartest2(categoryScores.spatialLag{1},categoryScores.spatialLag{2});
[h,p] = vartest2(categoryScores.randomMap{1},categoryScores.randomMap{2});

%-------------------------------------------------------------------------------
categoryScoresTogether = [categoryScores.randomMap; categoryScores.spatialLag];
categoryLabelsTogether = [categoryLabels.randomMap; categoryLabels.spatialLag];

%-------------------------------------------------------------------------------
% Violin plots:
f = figure('color','w');
extraParams = struct();
extraParams.theColors = {theColors(1,:),theColors(2,:),theColors(1,:),theColors(2,:)};
extraParams.customSpot = '';
extraParams.offsetRange = 0.7;
BF_JitteredParallelScatter(categoryScoresTogether,true,true,false,extraParams);
ax = gca();
ax.XTick = 1:4;
ax.XTickLabel = categoryLabelsTogether;
ylabel('Mean category score')
xlabel('GO category')
title(sprintf('%u nulls',numNullSamples))
f.Position = [1000        1078         341         260];
maxDev = max([abs(min([categoryScoresTogether{:}])),max([categoryScoresTogether{:}])]);
ax.YLim = [-maxDev,maxDev];
