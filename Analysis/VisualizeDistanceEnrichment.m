function VisualizeDistanceEnrichment(ResultsTable,categoryIndex,params)
% Visualize a specific GO category with respect to its spatial autocorrelation
% Visualize individual gene categories with different types of scores:
%-------------------------------------------------------------------------------

if nargin < 2
    categoryIndex = 1;
    % ResultsTable ordered—-higher means lower p-value
end
if nargin < 3
    params = GiveMeDefaultParams('mouse');
end

%-------------------------------------------------------------------------------
% Load in additional data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
distMat = GiveMeDistanceMatrix(params.humanOrMouse);

isThisCat = ismember(geneInfo.entrez_id,ResultsTable.annotations{categoryIndex});
geneDataSub = geneData(:,isThisCat);

%-------------------------------------------------------------------------------
f = figure('color','w');
f.Position = [619,679,1096,368];

% Clustered Data Matrix
ax = subplot(1,3,1);
ord_row = BF_ClusterReorder(geneDataSub,'corr','average');
ord_col = BF_ClusterReorder(geneDataSub','corr','average');
geneAcronyms = geneInfo.acronym(isThisCat);
geneAcronyms = geneAcronyms(ord_col);
BF_imagesc(BF_NormalizeMatrix(geneDataSub(ord_row,ord_col),'mixedSigmoid'))
colormap([flipud(BF_getcmap('blues',9));1,1,1;BF_getcmap('reds',9)])
ax.XTick = 1:size(geneDataSub,2);
ax.XTickLabel = geneAcronyms;
ax.XTickLabelRotation = 90;
ylabel('Genes')
ylabel('Brain areas')
title(ResultsTable.GOName{categoryIndex})

% Spatial projection:
subplot(1,3,2); axis('equal')
% nanmean(geneDataSub,2);
gMeanFlipped = flippedMean(geneDataSub);
gMean = BF_NormalizeMatrix(gMeanFlipped,'mixedSigmoid');
coOrdsXY = mdscale(distMat,2);
scatter(coOrdsXY(:,1),coOrdsXY(:,2),25,gMean,'filled')
colormap([flipud(BF_getcmap('blues',9));1,1,1;BF_getcmap('reds',9)])
xlabel('spatialAxis1')
ylabel('spatialAxis2')
cB = colorbar;
cB.Label.String = '(flipped) mean expression levels';

% Histogram of individual gene correlations:
ax = subplot(1,3,3); hold on
isUpperDiag = triu(true(size(distMat)),+1);
individualRhos = zeros(size(geneDataSub,2),1);
for j = 1:size(geneDataSub,2)
    GCCj = geneDataSub(:,j)*geneDataSub(:,j)';
    individualRhos(j) = corr(distMat(isUpperDiag),GCCj(isUpperDiag),'type','Spearman');
end
histogram(individualRhos)
GCCmean = gMean*gMean';
rhoMean = corr(distMat(isUpperDiag),GCCmean(isUpperDiag),'type','Spearman');
h1 = plot(ones(2,1)*rhoMean,ax.YLim,'r');
h2 = plot(-ones(2,1)*ResultsTable.meanScore(categoryIndex),ax.YLim,'g');
legend([h1,h2],'meanAgglomFlipped','meanIndividual')
xlabel('\rho (GCC-distance)')
% plot(distMat(isUpperDiag),GCC(isUpperDiag),'.k')
% ylabel('Mean expression GCC')
% title(sprintf('rho = %.3f',rho))

end
