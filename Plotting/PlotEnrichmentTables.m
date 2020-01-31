function [allGOIDsSort,GONamesSort] = PlotEnrichmentTables(resultsTables,thresholds,whatSpecies)
% Plots multiple results of GO enrichment as a (sorted) table
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check Inputs:
%-------------------------------------------------------------------------------
if nargin < 2
    % Need at least this many annotations to be included:
    thresholds = zeros(1,3);
    thresholds(1) = 0.2; % threshold for plotting a 'significant' category
    thresholds(2) = 3; % GO category needs at least 3 significant studies to be included
    thresholds(3) = 2; % Study needs to flag at least 2 significant categories to be included
end
sigThreshold = thresholds(1);
annotationThresholdGO = thresholds(2);
annotationThresholdDataset = thresholds(3);
maxCategories = 200; % show at most this many categories

if nargin < 3
    whatSpecies = '';
end

%-------------------------------------------------------------------------------
[rowVectorResults,allGOIDs,allTableNames] = CombineTables(resultsTables,whatSpecies,'pValZCorr');
numTables = length(allTableNames);

% Retrieve full GO information:
params = GiveMeDefaultParams(whatSpecies);
params.e.sizeFilter = [0,1e6];
GOTerms = GiveMeGOData(params);


% Set p-value for results over threshold to 1, in rowVectorResultsTh
underThreshold = @(x) 1-(1-x)*double(x < sigThreshold);
rowVectorResultsTh = arrayfun(underThreshold,rowVectorResults);
rowVectorResultsTh(isnan(rowVectorResultsTh)) = 1;

%-------------------------------------------------------------------------------
% TRIMMING
%-------------------------------------------------------------------------------
% rowVectorResultsTh has dimensions: [study x GO category]

% Remove studies with fewer than Y annotations:
numAnnotations = sum(rowVectorResultsTh < 1,2);
hasNoAnnotations = (numAnnotations < annotationThresholdDataset);
rowVectorResultsTh(hasNoAnnotations,:) = [];
allTableNames(hasNoAnnotations) = [];
numTrimmedDatasets = sum(~hasNoAnnotations);
fprintf(1,'Trimmed %u -> %u studies flagging fewer than %u annotations\n',...
        length(hasNoAnnotations),numTrimmedDatasets,annotationThresholdDataset);


% Trim out GO categories with fewer than X annotations:
numAnnotated = sum(rowVectorResultsTh < 1,1);
hasNoAnnotations = (numAnnotated < annotationThresholdGO);
allGOIDs(hasNoAnnotations) = [];
rowVectorResultsTh(:,hasNoAnnotations) = [];
numTrimmed = sum(~hasNoAnnotations);
fprintf(1,'Trimmed %u -> %u GO Categories with fewer than %u annotations\n',...
                    length(hasNoAnnotations),numTrimmed,annotationThresholdGO);

% Trim GO categories further based on a fixed max number:
if size(rowVectorResultsTh,2) > maxCategories
    fprintf(1,'Trimmed %u -> %u GO Categories by fixed limit\n',...
                                size(rowVectorResultsTh,2),maxCategories);
    meansGO = nanmean(rowVectorResultsTh,1);
    [~,ix] = sort(meansGO,'ascend');
    keepMe = ix(1:maxCategories);
    rowVectorResultsTh = rowVectorResultsTh(:,keepMe);
    allGOIDs = allGOIDs(keepMe);
    numTrimmed = maxCategories;
end

% Summarize each dataset and each GO category by mean enrichment:
meansGO = nanmean(rowVectorResultsTh,1);
meansDatasets = nanmean(rowVectorResultsTh,2);

%-------------------------------------------------------------------------------
% Sort rows/columns for visualization:
%-------------------------------------------------------------------------------
[~,ix] = sort(meansGO,'ascend');
allGOIDsSort = allGOIDs(ix);
% Sort rows by number of annotations, or by a linkage clustering:
% iy = BF_ClusterReorder(rowVectorResultsTh,'euclidean','average');
[~,iy] = sort(meansDatasets,'ascend');
allTableNamesSort = allTableNames(iy);

rowVectorResultsThSort = rowVectorResultsTh(iy,ix);

% IDs -> GONames
GONamesSort = cell(numTrimmed,1);
for i = 1:numTrimmed
    weHere = (GOTerms.GOID==allGOIDsSort(i));
    if ~any(weHere)
        GONamesSort{i} = sprintf('<unknown> (ID:%u)',allGOIDsSort(i));
    else
        GONamesSort{i} = sprintf('[%u]%s(%u)',GOTerms.GOID(weHere),GOTerms.GOName{weHere},GOTerms.size(weHere));
    end
    fprintf(1,'%u %s(%u)\n',allGOIDsSort(i),GOTerms.GOName{weHere},GOTerms.size(weHere));
end

%-------------------------------------------------------------------------------
% Plot:
f = figure('color','w'); ax = gca;
imagesc(rowVectorResultsThSort');
ax.XTick = 1:numTables;
ax.XTickLabel = allTableNamesSort;
ax.XTickLabelRotation = 40;
ax.TickLabelInterpreter = 'none';
ax.YTick = 1:sum(~hasNoAnnotations);
% ax.XTickLabel = allGOIDsSort;
ax.YTickLabel = GONamesSort;
colormap([flipud(BF_getcmap('purplebluegreen',9));1,1,1]);
caxis([0,sigThreshold*1.2]);
title(whatSpecies);

cB = colorbar();
cB.Label.String = 'pCorr';

end
