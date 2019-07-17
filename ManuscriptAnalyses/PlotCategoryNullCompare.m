% PlotCategoryNullCompare
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Find two good candidates:

% Load VE1 scores:
whatSpecies = 'mouse';
whatShuffle = 'geneShuffle'; % 'geneShuffle', 'independentSpatialShuffle'
numNullSamples_VE1 = 20000; % (Intra_*_VE1_20000.mat)
resultsIntra = load(sprintf('Intra_%s_%s_VE1_%u.mat',whatSpecies,whatShuffle,numNullSamples_VE1));
categorySizeRange = [40,42];
isInSizeRange = (resultsIntra.resultsTable.size >= categorySizeRange(1)) & (resultsIntra.resultsTable.size <= categorySizeRange(2))
resultsIntraSized = resultsIntra.resultsTable(isInSizeRange,:)

%-------------------------------------------------------------------------------
% Parameters:
% whatGOIDs = [7215,32612];
whatGOIDs = [31638,61001];
numNullSamples = 20000;
whatCorr = 'Spearman';

% Assign colors:
theColors = [[32,178,170]/255;[184,134,11]/255];
% [220,220,220]/255}; % [119,136,153]/255, [119,136,153]/255, [220,220,220]/255};


%-------------------------------------------------------------------------------
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

%-------------------------------------------------------------------------------
categoryScores = struct();
categoryLabels = struct();
fprintf(1,'Computing distribution of null spatial-lag category scores for %u GO categories\n',length(whatGOIDs));
[categoryScores.spatialLag,categoryLabels.spatialLag] = CompareNulls(whatGOIDs,whatSpecies,'spatialLag',numNullSamples,false);
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
extraParams.theColors = {theColors(1,:),theColors(2,:),brighten(theColors(1,:),+0.5),brighten(theColors(2,:),+0.5)};
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
