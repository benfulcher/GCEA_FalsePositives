function PlotExpressionMatrix(whatSpecies)
% Plot an expression matrix of a set of genes across brain areas
%-------------------------------------------------------------------------------

if nargin < 1
    whatSpecies = 'mouse';
end

params = GiveMeDefaultParams(whatSpecies);
params.g.normalizationGene = 'scaledRobustSigmoid';
params.g.normalizationRegion = 'none';
[geneData,geneInfo,structInfo] = LoadMeG(params.g);

%-------------------------------------------------------------------------------
% Cluster-reorder genes
ord_col = BF_ClusterReorder(geneData','euclidean','average');

%-------------------------------------------------------------------------------
% Plot gene-reordered matrix:
PlotColorMatrix(geneData(:,ord_col),structInfo);

end
