% function LiteratureLook(whatSpecies)

% if nargin < 1
    whatSpecies = 'mouse';
% end

% Load in data:
load('LiteratureEnrichmentLoaded.mat','resultsTables','mouseOrHuman');

% Filter by species:
allTableNames = fieldnames(resultsTables);
theSpeciesTables = find(mouseOrHuman==whatSpecies);
speciesTableNames = allTableNames(theSpeciesTables);
numSpecies = length(theSpeciesTables);
resultsTablesSpecies = struct;
for i = 1:numSpecies
    theName = speciesTableNames{i};
    resultsTablesSpecies.(theName) = resultsTables.(theName);
end
fprintf(1,'%u GO enrichment datasets involving %s\n',numSpecies,whatSpecies);

% Get all GOIDs
allGOIDs = structfun(@(x)x.GOID,resultsTablesSpecies,'UniformOutput',false);
allGOIDs = struct2cell(allGOIDs);
allGOIDs = sort(unique(vertcat(allGOIDs{:})),'ascend');
numGOIDs = length(allGOIDs);
fprintf(1,'Annotations span %u GO categories\n',length(allGOIDs));

% Convert each dataset to a vector across GOIDs:
rowVectorResults = nan(numSpecies,numGOIDs);
for i = 1:numSpecies
    [~,ia,ib] = intersect(allGOIDs,resultsTablesSpecies.(speciesTableNames{i}).GOID);
    rowVectorResults(i,ia) = resultsTablesSpecies.(speciesTableNames{i}).pValCorr(ib);
end

% Filter max:
theThreshold = 0.1;
underThreshold = @(x) 1-(1-x)*double(x < theThreshold); % results over threshold are set to 1
rowVectorResultsTh = arrayfun(underThreshold,rowVectorResults);
rowVectorResultsTh(rowVectorResultsTh==1) = NaN;

% Sort:
meansGO = nanmean(~isnan(rowVectorResultsTh),1);
meansDatasets = nanmean(~isnan(rowVectorResultsTh),2);

[~,ix] = sort(meansGO,'descend');
allGOIDsSort = allGOIDs(ix);
[~,iy] = sort(meansDatasets,'descend');
speciesTableNamesSort = speciesTableNames(iy);

rowVectorResultsThSort = rowVectorResultsTh(iy,ix);

% Trim:
hasNoAnnotations = (meansGO(ix)==0);
allGOIDsSort(hasNoAnnotations) = [];
rowVectorResultsThSort(:,hasNoAnnotations) = [];

% Plot:
f = figure('color','w'); ax = gca;
BF_imagesc(rowVectorResultsThSort);
ax.YTick = 1:numSpecies;
ax.YTickLabel = flipud(speciesTableNamesSort);
ax.XTick = 1:sum(~hasNoAnnotations);
ax.XTickLabel = allGOIDsSort;
colormap([flipud(BF_getcmap('blues',9))])
title(whatSpecies);

% end
