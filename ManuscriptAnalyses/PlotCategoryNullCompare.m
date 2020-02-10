% function PlotCategoryNullCompare(whatSpecies)
% if nargin < 1
whatSpecies = 'mouse';
% end
%-------------------------------------------------------------------------------

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
% ---Plot information about the specific categories (and a ranodm one)---
doReorder = false;
newFigure = false;
doLabels = false;
categoryWhat = 1;
f = figure('color','w');
%-------------------------------------------------------------------------------
% CATEGORY 1:
subplot(5,3,[1,4,7])
PlotCategoryExpression(whatGOIDs(1),params,doReorder,newFigure,doLabels);
subplot(5,3,[10,13])
PlotCategoryIntraCorr(whatGOIDs(1),params,params.e.whatCorr,newFigure,doLabels);
%-------------------------------------------------------------------------------
% CATEGORY 2:
subplot(5,3,[2,5,8])
PlotCategoryExpression(whatGOIDs(2),params,doReorder,newFigure,doLabels);
subplot(5,3,[11,14])
PlotCategoryIntraCorr(whatGOIDs(2),params,params.e.whatCorr,newFigure,doLabels);
%-------------------------------------------------------------------------------
% CATEGORY 3:
subplot(5,3,[3,6,9])
PlotCategoryExpression([],params,doReorder,newFigure,doLabels);
subplot(5,3,[12,15])
PlotCategoryIntraCorr([],params,params.e.whatCorr,newFigure,doLabels);
f.Position = [1111         542         560         413];
%-------------------------------------------------------------------------------
% Save:
fileName = fullfile('OutputPlots','ExpressionCoexpression.svg');
saveas(f,fileName,'svg')
fprintf(1,'Saved to %s\n',fileName);

%===============================================================================

% -------------------------------------------------------------------------------
% FPSR of these categories?
FPSR_random = SurrogateEnrichmentProcess(whatSpecies,params.nulls.numNullsFPSR,'randomUniform','');
FPSR_spatial = SurrogateEnrichmentProcess(whatSpecies,params.nulls.numNullsFPSR,'spatialLag','');
itsMeRand = find(ismember(FPSR_random.GOID,whatGOIDs));
itsMeSpat = find(ismember(FPSR_spatial.GOID,whatGOIDs));
for i = 1:length(itsMeRand)
    fprintf(1,'Category ''%s'' has SBP-random FPSR of %.3f%%\n',...
                            FPSR_random.GOName{itsMeRand(i)},...
                            FPSR_random.sumUnderSig(itsMeRand(i))/1e4*100);
    fprintf(1,'Category ''%s'' has SBP-spatial FPSR of %.3f%%\n',...
                            FPSR_spatial.GOName{itsMeSpat(i)},...
                            FPSR_spatial.sumUnderSig(itsMeSpat(i))/1e4*100);
end
%===============================================================================

%-------------------------------------------------------------------------------
% Now let's get going with a joint null samples for this category size:
%-------------------------------------------------------------------------------
numNullSamplesJoint = params.e.numNullSamples;
params.e.whatEnsemble = 'randomMap';
nullDistributionRandomMap = ComputeJointNull(params,40,numNullSamplesJoint);
params.e.whatEnsemble = 'customEnsemble';
nullDistributionSpatialLag = ComputeJointNull(params,40,numNullSamplesJoint);
% f = figure('color','w');
% hold('on')
% histogram(nullDistributionRandomMap)
% histogram(nullDistributionSpatialLag)

%===============================================================================
% Plot FPSR distributions
%-------------------------------------------------------------------------------
categoryScores = struct();
categoryLabels = struct();
[categoryScores.spatialLag,categoryLabels.spatialLag] = CompareNulls(whatGOIDs,whatSpecies,'customEnsemble',false);
[categoryScores.randomMap,categoryLabels.randomMap] = CompareNulls(whatGOIDs,whatSpecies,'randomMap',false);

%-------------------------------------------------------------------------------
% Tests of variance just for fun:
[h,p] = vartest2(categoryScores.spatialLag{1},categoryScores.spatialLag{2});
[h,p] = vartest2(categoryScores.randomMap{1},categoryScores.randomMap{2});

%-------------------------------------------------------------------------------
% Add random-gene joint nulls?:
categoryScores.spatialLag{3} = nullDistributionSpatialLag;
categoryScores.randomMap{3} = nullDistributionRandomMap;
categoryLabels.spatialLag{3} = 'Random Gene';
categoryLabels.randomMap{3} = 'Random Gene';

%-------------------------------------------------------------------------------
categoryScoresTogether = [categoryScores.randomMap; categoryScores.spatialLag];
categoryLabelsTogether = [categoryLabels.randomMap; categoryLabels.spatialLag];

%-------------------------------------------------------------------------------
% Violin plots:
theColors = GiveMeColors('twoGOCategories'); % (assign colors)
f = figure('color','w');
extraParams = struct();
extraParams.theColors = {theColors(1,:),theColors(2,:),theColors(3,:),...
                            theColors(1,:),theColors(2,:),theColors(3,:)};
extraParams.customSpot = '';
extraParams.offsetRange = 0.7;
BF_JitteredParallelScatter(categoryScoresTogether,true,true,false,extraParams);
ax = gca();
ax.XTick = 1:6;
ax.XTickLabel = categoryLabelsTogether;
ylabel('Mean category score')
xlabel('GO category')
title(sprintf('%u nulls',params.e.numNullSamples))
f.Position = [1000        1078         341         260];
maxDev = max([abs(min([categoryScoresTogether{:}])),max([categoryScoresTogether{:}])]);
ax.YLim = [-maxDev,maxDev];

% Save out as svg file:
fileName = fullfile('OutputPlots','NullCompareDistributions.svg');
saveas(f,fileName,'svg')
fprintf(1,'Saved to %s\n',fileName);
