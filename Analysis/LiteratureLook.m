% function LiteratureLook(whatSpecies)

% if nargin < 1
    whatSpecies = 'human';
    theThreshold = 0.2;
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

% Match to names
eParam = GiveMeDefaultParams('enrichment');
eParam.sizeFilter = [0,1e6];
[GOTerms,geneEntrezAnnotations] = GetFilteredGOData(eParam.whatSource,...
                    eParam.processFilter,eParam.sizeFilter);

% Convert each dataset to a vector across GOIDs:
rowVectorResults = nan(numSpecies,numGOIDs);
for i = 1:numSpecies
    [~,ia,ib] = intersect(allGOIDs,resultsTablesSpecies.(speciesTableNames{i}).GOID);
    rowVectorResults(i,ia) = resultsTablesSpecies.(speciesTableNames{i}).pValCorr(ib);
end

% Filter max:
underThreshold = @(x) 1-(1-x)*double(x < theThreshold); % results over threshold are set to 1
rowVectorResultsTh = arrayfun(underThreshold,rowVectorResults);
rowVectorResultsTh(isnan(rowVectorResultsTh)) = 1;

% Sort:
meansGO = nanmean(rowVectorResultsTh,1);
meansDatasets = nanmean(rowVectorResultsTh,2);

[~,ix] = sort(meansGO,'ascend');
allGOIDsSort = allGOIDs(ix);
[~,iy] = sort(meansDatasets,'ascend');
speciesTableNamesSort = speciesTableNames(iy);

rowVectorResultsThSort = rowVectorResultsTh(iy,ix);

% Trim:
numTrimmed = sum(~hasNoAnnotations);
hasNoAnnotations = (meansGO(ix)==1);
allGOIDsSort(hasNoAnnotations) = [];
rowVectorResultsThSort(:,hasNoAnnotations) = [];

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

% Plot:
f = figure('color','w'); ax = gca;
imagesc(rowVectorResultsThSort);
ax.YTick = 1:numSpecies;
ax.YTickLabel = flipud(speciesTableNamesSort);
ax.TickLabelInterpreter = 'none';
ax.XTick = 1:sum(~hasNoAnnotations);
% ax.XTickLabel = allGOIDsSort;
ax.XTickLabel = GONamesSort;
ax.XTickLabelRotation = 90;
colormap([flipud(BF_getcmap('blues',9))])
caxis([0,theThreshold*1.2]);
title(whatSpecies);

% end
