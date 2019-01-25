function PlotExpressionMatrix(whatSpecies)
% Plot an expression matrix of a set of genes
%-------------------------------------------------------------------------------

if nargin < 1
    whatSpecies = 'mouse';
end

params = GiveMeDefaultParams(whatSpecies);
params.g.normalizationGene = 'scaledSigmoid';
[geneData,geneInfo,structInfo] = LoadMeG(params.g);

%-------------------------------------------------------------------------------
% Cluster-reorder genes
ord_col = BF_ClusterReorder(geneData','corr','average');

%-------------------------------------------------------------------------------
% Plot:
PlotColorMatrix(geneData(:,ord_col),structInfo);

end
