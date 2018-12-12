function PlotCategoryExpression(whatGOID,params)
% Plots the spatial expression patterns of genes within a given GO category
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs:
if nargin < 1
    whatGOID = 6099;
end
if nargin < 2
    params = GiveMeDefaultParams('mouse');
end

%-------------------------------------------------------------------------------
% Load in gene expression data for this category
[geneData,geneInfo,structInfo,categoryInfo] = GiveMeGOCategory(whatGOID,params);
numGenes = height(geneInfo);

%-------------------------------------------------------------------------------
% Plot a clustered gene-expression data matrix
f = figure('color','w');
f.Position = [619,679,1096,368];
ax = gca;
ord_row = BF_ClusterReorder(geneData,'corr','average');
ord_col = BF_ClusterReorder(geneData','corr','average');
BF_imagesc(BF_NormalizeMatrix(geneData(ord_row,ord_col),'mixedSigmoid'))
colormap([flipud(BF_getcmap('blues',9));1,1,1;BF_getcmap('reds',9)])
ax.XTick = 1:numGenes;
ax.XTickLabel = geneInfo.acronym(ord_col);
ax.XTickLabelRotation = 90;
ylabel('Genes')
ylabel('Brain areas')
title(categoryInfo.GOName{1})

end
