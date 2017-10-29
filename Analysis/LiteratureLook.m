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

% end
