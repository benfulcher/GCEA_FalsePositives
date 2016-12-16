function [GOTable,geneEntrezAnnotations] = GetFilteredGOData(whatFilter,sizeFilter)

if nargin < 1
    whatFilter = 'biological_process';
end
if nargin < 2
    sizeFilter = [5,200];
end

% Get GO annotation data (processed):
load('GOAnnotation.mat','allGOCategories','geneEntrezAnnotations');

% Get GO ontology details
GOTable = GetGOTerms(whatFilter);

%-------------------------------------------------------------------------------
% Filter
%-------------------------------------------------------------------------------
% Filter by ontology details:
[~,ia,ib] = intersect(GOTable.GOID,allGOCategories);
fprintf(1,'Filtering to %u annotated GO categories related to %s\n',length(ia),whatFilter);
GOTable = GOTable(ia,:);
allGOCategories = allGOCategories(ib);
geneEntrezAnnotations = geneEntrezAnnotations(ib);

% Filter by category size:
numGOCategories = length(allGOCategories);
sizeGOCategories = cellfun(@length,geneEntrezAnnotations);
isGoodSize = (sizeGOCategories >= sizeFilter(1)) & (sizeGOCategories <= sizeFilter(2));
allGOCategories = allGOCategories(isGoodSize);
geneEntrezAnnotations = geneEntrezAnnotations(isGoodSize);
GOTable = GOTable(isGoodSize,:);
sizeGOCategories = sizeGOCategories(isGoodSize);
numGOCategories = length(allGOCategories);
fprintf(1,'Filtered to %u categories with between %u and %u annotations\n',...
                numGOCategories,sizeFilter(1),sizeFilter(2));

end
