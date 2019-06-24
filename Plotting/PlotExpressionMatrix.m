function PlotExpressionMatrix(whatSpecies,doRandomize)
% Plot an expression matrix of a set of genes across brain areas
%-------------------------------------------------------------------------------

if nargin < 1
    whatSpecies = 'mouse';
end
if nargin < 2
    doRandomize = false;
end

%-------------------------------------------------------------------------------
params = GiveMeDefaultParams(whatSpecies);
params.g.normalizationGene = 'scaledRobustSigmoid';
params.g.normalizationRegion = 'none';
[geneData,geneInfo,structInfo] = LoadMeG(params.g);

%-------------------------------------------------------------------------------
% Cluster-reorder genes
if doRandomize
    numGenes = size(geneData,2);
    numAreas = size(geneData,1);
    for i = 1:numGenes
        rp = randperm(numAreas);
        geneData(:,i) = geneData(rp,i);
    end
end

ord_col = BF_ClusterReorder(geneData','euclidean','average');

%-------------------------------------------------------------------------------
% Plot gene-reordered matrix:
PlotColorMatrix(geneData(:,ord_col),structInfo);

end
