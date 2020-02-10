function PlotCategoryExpression(whatGOID,params,doReorder,newFigure,doLabels)
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
if nargin < 3
    doReorder = true;
end
if nargin < 4
    newFigure = true;
end
if nargin < 5
    doLabels = true;
end

%-------------------------------------------------------------------------------
% Load in gene expression data for this category
[geneData,geneInfo,structInfo,categoryInfo] = GiveMeGOCategory(whatGOID,params);
numGenes = height(geneInfo);

%-------------------------------------------------------------------------------
% Plot a clustered gene-expression data matrix
if newFigure
    f = figure('color','w');
    f.Position = [619,679,1096,368];
end
ax = gca;
ord_col = BF_ClusterReorder(geneData','corr','average');
if doReorder
    ord_row = BF_ClusterReorder(geneData,'corr','average');
else
    ord_row = 1:size(geneData,1);
end
BF_imagesc(BF_NormalizeMatrix(geneData(ord_row,ord_col),'mixedSigmoid'))
colormap([flipud(BF_getcmap('blues',9));1,1,1;BF_getcmap('reds',9)])
ax.XTick = 1:numGenes;
if doLabels
    ax.XTickLabel = geneInfo.acronym(ord_col);
    ax.XTickLabelRotation = 90;
    xlabel('Genes')
    ylabel('Brain areas')
    try
        title(categoryInfo.GOName{1})
    end
else
    ax.XTickLabel = '';
    ax.YTickLabel = '';
end

end
