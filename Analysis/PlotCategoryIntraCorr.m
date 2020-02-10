function PlotCategoryIntraCorr(whatGOID,params,whatCorr,newFigure,doLabeling)
% Plots the intra-category correlations between genes
%-------------------------------------------------------------------------------

% Check inputs:
if nargin < 1
    whatGOID = 6099;
end
if nargin < 2
    params = GiveMeDefaultParams('mouse');
end
if nargin < 3
    whatCorr = 'Spearman';
end
if nargin < 4
    newFigure = true;
end
if nargin < 5
    doLabeling = true;
end

%-------------------------------------------------------------------------------
% Load in gene expression data for this category
[geneData,geneInfo,structInfo,categoryInfo] = GiveMeGOCategory(whatGOID,params);
numGenes = height(geneInfo);

%-------------------------------------------------------------------------------
if newFigure
    f = figure('color','w');
end
ax = gca;
corrMat = corr(geneData,'type',whatCorr,'rows','pairwise');
ord = BF_ClusterReorder(corrMat,'euclidean','average');
imagesc(corrMat(ord,ord));

if doLabeling
    ax.YTick = 1:numGenes;
    ax.YTickLabel = geneInfo.acronym(ord);
    try
        title(sprintf('%s %s (%u genes)[%s]',categoryInfo.GOIDlabel{1},categoryInfo.GOName{1},....
                                    categoryInfo.size,whatCorr))
    end
    cB = colorbar;
else
    ax.YTickLabel = '';
    ax.XTickLabel = '';
end

caxis([-1,1])
colormap([flipud(BF_getcmap('blues',9)); 1,1,1; BF_getcmap('reds',9)])
axis('square')

end
