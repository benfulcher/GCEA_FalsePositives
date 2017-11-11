function PlotEnrichmentTables(resultsTables,theThreshold,whatSpecies)
% Plots multiple results of GO enrichment as a (sorted) table
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check Inputs:
%-------------------------------------------------------------------------------
if nargin < 2
    theThreshold = 0.2; % threshold for plotting a 'significant' category
end
if nargin < 3
    whatSpecies = '';
end
maxCategories = 200; % show at most this many categories
%-------------------------------------------------------------------------------

% Preliminaries:
allTableNames = fieldnames(resultsTables);
numTables = length(allTableNames);
fprintf(1,'%u GO enrichment datasets\n',numTables);

% Get all GOIDs
allGOIDs = structfun(@(x)x.GOID,resultsTables,'UniformOutput',false);
allGOIDs = struct2cell(allGOIDs);
allGOIDs = sort(unique(vertcat(allGOIDs{:})),'ascend');
numGOIDs = length(allGOIDs);
fprintf(1,'Annotations span %u GO categories\n',length(allGOIDs));

% Match to names
params = GiveMeDefaultParams(whatSpecies);
eParam = params.e;
eParam.sizeFilter = [0,1e6];
GOTerms = GetFilteredGOData(eParam.whatSource,eParam.processFilter,eParam.sizeFilter);

% Convert each dataset to a vector across GOIDs:
rowVectorResults = nan(numTables,numGOIDs);
for i = 1:numTables
    [~,ia,ib] = intersect(allGOIDs,resultsTables.(allTableNames{i}).GOID);
    rowVectorResults(i,ia) = resultsTables.(allTableNames{i}).pValCorr(ib);
end

% Filter max:
underThreshold = @(x) 1-(1-x)*double(x < theThreshold); % results over threshold are set to 1
rowVectorResultsTh = arrayfun(underThreshold,rowVectorResults);
rowVectorResultsTh(isnan(rowVectorResultsTh)) = 1;

% Summarize each dataset and each GO category by mean enrichment:
meansGO = nanmean(rowVectorResultsTh,1);
meansDatasets = nanmean(rowVectorResultsTh,2);

%-------------------------------------------------------------------------------
% TRIMMING
%-------------------------------------------------------------------------------
% Need at least this many annotations to be included:
annotationThresholdDataset = 3; % at least 3 significant to be included
annotationThresholdGO = 2; % category implicated in at least two studies to be included

% Trim out studies with fewer than Y annotations:
numAnnotations = sum(rowVectorResultsTh < 1,2);
hasNoAnnotations = (numAnnotations < annotationThresholdDataset);
allTableNames(hasNoAnnotations) = [];
rowVectorResultsTh(hasNoAnnotations,:) = [];

numTrimmedDatasets = sum(~hasNoAnnotations);
fprintf(1,'Trimmed %u -> %u datasets with fewer than %u annotations\n',...
            length(hasNoAnnotations),numTrimmedDatasets,annotationThresholdDataset);

% Trim out GO categories with fewer than X annotations:
numAnnotated = sum(rowVectorResultsTh < 1,1);
hasNoAnnotations = (numAnnotated < annotationThresholdGO);
allGOIDs(hasNoAnnotations) = [];
rowVectorResultsTh(:,hasNoAnnotations) = [];
meansGO(hasNoAnnotations) = [];
numTrimmed = sum(~hasNoAnnotations);
fprintf(1,'Trimmed %u -> %u GO Categories with fewer than %u annotations\n',...
                    length(hasNoAnnotations),numTrimmed,annotationThresholdGO);

% Trim further based on the number
if size(rowVectorResultsTh,2) > maxCategories
    rowVectorResultsTh = rowVectorResultsTh(:,1:maxCategories);
    meansGO = meansGO(1:maxCategories);
    numTrimmed = maxCategories;
end

%-------------------------------------------------------------------------------
% Sort rows/columns for visualization:
%-------------------------------------------------------------------------------
[~,ix] = sort(meansGO,'ascend');
allGOIDsSort = allGOIDs(ix);
% Sort rows by number of annotations, or by a linkage clustering:
iy = BF_ClusterReorder(rowVectorResultsTh,'euclidean','average');
% [~,iy] = sort(meansDatasets,'ascend');
allTableNamesSort = allTableNames(iy);

rowVectorResultsThSort = rowVectorResultsTh(iy,ix);

% IDs -> GONames
GONamesSort = cell(numTrimmed,1);
for i = 1:numTrimmed
    weHere = (GOTerms.GOID==allGOIDsSort(i));
    if ~any(weHere)
        GONamesSort{i} = '<unknown>';
    else
        GONamesSort{i} = sprintf('%s(%u)',GOTerms.GOName{weHere},GOTerms.size(weHere));
    end
end

%-------------------------------------------------------------------------------
% Plot:
f = figure('color','w'); ax = gca;
imagesc(rowVectorResultsThSort);
ax.YTick = 1:numTables;
ax.YTickLabel = allTableNamesSort;
ax.TickLabelInterpreter = 'none';
ax.XTick = 1:sum(~hasNoAnnotations);
% ax.XTickLabel = allGOIDsSort;
ax.XTickLabel = GONamesSort;
ax.XTickLabelRotation = 90;
colormap([flipud(BF_getcmap('blues',9))])
caxis([0,theThreshold*1.2]);
title(whatSpecies)

end
