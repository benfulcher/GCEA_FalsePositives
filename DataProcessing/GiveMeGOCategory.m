function [geneDataSub,geneInfoSub,structInfo,categoryInfo] = GiveMeGOCategory(whatGOID,params)
% Retrieve expression data for a given GO category
if nargin < 1
    whatGOID = 7215;
end
if nargin < 2
    params = GiveMeDefaultParams('mouse');
end

%-------------------------------------------------------------------------------
% Load gene-expression data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
% Load GO annotations:
GOTable = GetFilteredGOData(params.e.dataSource,params.e.processFilter,...
                                params.e.sizeFilter,geneInfo.entrez_id);
% Get the category of interest:
whatCategory = find(GOTable.GOID==whatGOID);
% Match genes:
theGenesEntrez = GOTable.annotations{whatCategory};
[entrezMatched,ia,ib] = intersect(theGenesEntrez,geneInfo.entrez_id);
fprintf(1,'Looking in at %s:%s (%u)\n',GOTable.GOIDlabel{whatCategory},...
                    GOTable.GOName{whatCategory},GOTable.size(whatCategory));
fprintf(1,'%u/%u genes from this GO category match\n',length(entrezMatched),length(theGenesEntrez));
numGenesGO = length(entrezMatched);
matchMe = find(ismember(geneInfo.entrez_id,entrezMatched));

% Get set of genes corresponding to this category:
geneDataSub = geneData(:,matchMe);
geneInfoSub = geneInfo(matchMe,:);

% Information about the category:
categoryInfo = GOTable(whatCategory,:);

end
