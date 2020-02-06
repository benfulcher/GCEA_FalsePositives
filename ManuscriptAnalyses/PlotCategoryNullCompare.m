function PlotCategoryNullCompare(whatSpecies)
if nargin < 1
    whatSpecies = 'mouse';
end

% Set the GOIDs of the categories you want to look into:
% whatGOIDs = [7215,32612];
whatGOIDs = [61001,31638];

%-------------------------------------------------------------------------------
% Looks in detail at the null score distributions of two selected GO categories

params = GiveMeDefaultParams(whatSpecies);
% Also get some plots about one of the categories
plotInfoAboutSpecific = false;

%-------------------------------------------------------------------------------
% x x x
% | | | HOW THE TWO CANDIDATES WERE SELECTED | | |
% x x x
%-------------------------------------------------------------------------------
% Find two good candidate GO categories:
% Load intra-category correlation scores:
% whatSpecies = 'mouse';
% whatShuffle = 'geneShuffle'; % 'geneShuffle', 'independentSpatialShuffle'
% whatIntraStat = 'raw';
% numNullSamples_intraCorr = 20000; % (Intra_*_*_20000.mat)
% fileNameIntra = sprintf('Intra_%s_%s_%s_%u.mat',whatSpecies,whatShuffle,whatIntraStat,numNullSamples_intraCorr);
% resultsIntra = load(fileNameIntra);
% categorySizeRange = [40,42];
% isInSizeRange = (resultsIntra.resultsTable.size >= categorySizeRange(1)) & (resultsIntra.resultsTable.size <= categorySizeRange(2))
% resultsIntraSized = resultsIntra.resultsTable(isInSizeRange,:)


%===============================================================================
% ---Plot information about a specific category---
if plotInfoAboutSpecific
    categoryWhat = 2;
    PlotCategoryIntraCorr(whatGOIDs(categoryWhat),params,params.e.whatCorr);
    PlotCategoryExpression(whatGOIDs(categoryWhat),params);

    % -------------------------------------------------------------------------------
    % FPSR for this category?
    FPSR_random = SurrogateEnrichmentProcess(whatSpecies,params.nulls.numNullsFPSR,'randomUniform','');
    FPSR_thisCategory = FPSR_random.sumUnderSig(FPSR_random.GOID==whatGOIDs(categoryWhat));
    FPSR_random(4285,:)
end
%===============================================================================

%===============================================================================
% Plot FPSR distributions
%-------------------------------------------------------------------------------
categoryScores = struct();
categoryLabels = struct();
fprintf(1,'Computing distribution of null spatial-lag category scores for %u GO categories\n',...
                                                length(whatGOIDs));
[categoryScores.spatialLag,categoryLabels.spatialLag] = CompareNulls(whatGOIDs,whatSpecies,'customEnsemble',false);
fprintf(1,'Computing distribution of null random-map category scores for %u GO categories\n',...
                                                length(whatGOIDs));
[categoryScores.randomMap,categoryLabels.randomMap] = CompareNulls(whatGOIDs,whatSpecies,'randomMap',false);

%-------------------------------------------------------------------------------
% Tests of variance:
[h,p] = vartest2(categoryScores.spatialLag{1},categoryScores.spatialLag{2});
[h,p] = vartest2(categoryScores.randomMap{1},categoryScores.randomMap{2});

%-------------------------------------------------------------------------------
categoryScoresTogether = [categoryScores.randomMap; categoryScores.spatialLag];
categoryLabelsTogether = [categoryLabels.randomMap; categoryLabels.spatialLag];

%-------------------------------------------------------------------------------
% Violin plots:
theColors = GiveMeColors('twoGOCategories'); % (assign colors)
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
title(sprintf('%u nulls',params.e.numNullSamples))
f.Position = [1000        1078         341         260];
maxDev = max([abs(min([categoryScoresTogether{:}])),max([categoryScoresTogether{:}])]);
ax.YLim = [-maxDev,maxDev];
