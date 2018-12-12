%-------------------------------------------------------------------------------
% CompareNulls
%-------------------------------------------------------------------------------
% Plots comparison of the correlation between the genes in a given GO category
% and different null spatial maps
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Settings:
whatGOID = 6099;
whatSpecies = 'mouse';
numNulls = 100;

%-------------------------------------------------------------------------------
% Load in gene-expression data:
params = GiveMeDefaultParams(whatSpecies);
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
% Load in GO annotations:
GOTable = GetFilteredGOData(params.e.dataSource,params.e.processFilter,params.e.sizeFilter,...
                                    geneInfo.entrez_id);
numGOCategories = height(GOTable);
% Get the category of interest:
whatCategory = find(GOTable.GOID==whatGOID);
fprintf(1,'Looking at %s:%s (%u)\n',GOTable.GOIDlabel{whatCategory},...
                    GOTable.GOName{whatCategory},GOTable.size(whatCategory));

%-------------------------------------------------------------------------------
% Generate lots of null spatial maps:
