function [geneDataSub,geneInfoSub,structInfo,categoryInfo] = GiveMeGOCategory(whatGOID,params)
% Retrieve gene expression data for a given GO category
%-------------------------------------------------------------------------------
if nargin < 1
    whatGOID = 7215;
end
if nargin < 2
    params = GiveMeDefaultParams('mouse');
end

%-------------------------------------------------------------------------------
% Load gene-expression data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);

if isempty(whatGOID)
    warning('Returning a random set of 40 genes! That''s all!')
    numGenes = height(geneInfo);
    rp = randperm(numGenes);
    rpFilt = rp(1:40);
    geneDataSub = geneData(:,rpFilt);
    geneInfoSub = geneInfo(rpFilt,:);
    categoryInfo = table();
    return
end

% Load GO annotations to retrieve the category of interest:
GOTable = GiveMeGOData(params,geneInfo.entrez_id);
whatCategory = find(GOTable.GOID==whatGOID);
if isempty(whatCategory)
    error('No matches for category %u',whatGOID);
end

% Match genes:
theGenesEntrez = GOTable.annotations{whatCategory};
[entrezMatched,ia,ib] = intersect(theGenesEntrez,geneInfo.entrez_id);
fprintf(1,'Looking in at %s:%s (%u)\n',GOTable.GOIDlabel{whatCategory},...
                    GOTable.GOName{whatCategory},GOTable.size(whatCategory));
fprintf(1,'%u/%u genes from this GO category have matching records in the expression data\n',...
                        length(entrezMatched),length(theGenesEntrez));
numGenesGO = length(entrezMatched);
matchMe = find(ismember(geneInfo.entrez_id,entrezMatched));

% Get set of genes corresponding to this category:
geneDataSub = geneData(:,matchMe);
geneInfoSub = geneInfo(matchMe,:);

% Information about the category:
categoryInfo = GOTable(whatCategory,:);

end
