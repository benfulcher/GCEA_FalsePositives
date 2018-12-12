function PlotCategoryIntraCorr(whatGOID,params,whatCorr)
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

%-------------------------------------------------------------------------------
% Load in gene expression data for this category
[geneData,geneInfo,structInfo,categoryInfo] = GiveMeGOCategory(whatGOID,params);
numGenes = height(geneInfo);

%-------------------------------------------------------------------------------
f = figure('color','w');
ax = gca;
corrMat = corr(geneData,'type',whatCorr,'rows','pairwise');
ord = BF_ClusterReorder(corrMat,'euclidean','average');
imagesc(corrMat(ord,ord));
ax.YTick = 1:numGenes;
ax.YTickLabel = geneInfo.acronym(ord);
caxis([-1,1])
colormap([flipud(BF_getcmap('blues',9)); 1,1,1; BF_getcmap('reds',9)])
title(sprintf('%s %s (%u genes)[%s]',categoryInfo.GOIDlabel{1},categoryInfo.GOName{1},....
                                categoryInfo.size,whatCorr))
axis('square')
cB = colorbar;

end
