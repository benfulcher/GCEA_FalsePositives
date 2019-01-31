function [rowVectorResults,GOTerms,allGOIDs,allTableNames] = CombineTables(resultsTables,whatSpecies,customField)
% Structure of results tables (from GO enrichment analysis)
%-------------------------------------------------------------------------------

if nargin < 2
    whatSpecies = 'mouse';
end
if nargin < 3
    customField = 'pValCorr';
end
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
params.e.sizeFilter = [0,1e6];
GOTerms = GiveMeGOData(params);

% Convert each dataset to a vector across GOIDs:
rowVectorResults = nan(numTables,numGOIDs);
for i = 1:numTables
    [~,ia,ib] = intersect(allGOIDs,resultsTables.(allTableNames{i}).GOID);
    if iscell(customField)
        theField = customField{i};
    else
        theField = customField;
    end
    rowVectorResults(i,ia) = resultsTables.(allTableNames{i}).(theField)(ib);
end

end
