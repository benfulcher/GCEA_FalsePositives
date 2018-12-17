function [categoryScores,categoryLabels] = CompareNulls(whatGOIDs,whatSpecies,whatSurrogate,numNullSamples)
%-------------------------------------------------------------------------------
% Plots comparison of the correlation between the genes in a given GO category
% and different null spatial maps
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Settings:
if nargin < 1
    whatGOIDs = [7215,32612]; % 6099, 2374,2693,
end
if nargin < 2
    whatSpecies = 'mouse';
end
if nargin < 3
    % whatSurrogate = 'spatialShuffle';
    whatSurrogate = 'spatialLag';
end
if nargin < 4
    numNullSamples = 2000;
end
whatCorr = 'Spearman';
doPlot = false;

%-------------------------------------------------------------------------------
numGOIDs = length(whatGOIDs);
params = GiveMeDefaultParams(whatSpecies);
params.g.whatSurrogate = whatSurrogate;

%-------------------------------------------------------------------------------
% Loop over categories and compute null distributions:
categoryScores = cell(numGOIDs,1);
categoryLabels = cell(numGOIDs,1);
for i = 1:numGOIDs
    [categoryScores{i},categoryInfo] = GiveMeCategoryNullDist(whatGOIDs(i),...
                                                params,numNullSamples,whatCorr);
    categoryLabels{i} = categoryInfo.GOName{1};
end

%===============================================================================
if doPlot
    %-------------------------------------------------------------------------------
    % Violin plots:
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

    %-------------------------------------------------------------------------------
    % Plot null distributions as histograms:
    myColors = BF_getcmap('set2',6,0,0);

    f = figure('color','w');
    hold('on')
    h = cell(numGOIDs,1);
    for i = 1:numGOIDs
        h{i} = histogram(categoryScores{i},'normalization','pdf');
        h{i}.FaceColor = myColors(i,:);
    end
    legend([h{:}],categoryLabels)
    xlabel(sprintf('%s correlation',whatCorr))
    title(sprintf('%u nulls of %s',numNullSamples,whatSurrogate))
end
end
