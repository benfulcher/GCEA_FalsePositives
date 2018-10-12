% DistanceConfoundResults   Investigate the enrichment signatures of distance-related confounds

% Store results tables in this struct:
results = struct();

%-------------------------------------------------------------------------------
% Mouse:
params = GiveMeDefaultParams('mouse');
% To make GCC scores make sense -- expression needs to be [0,1] normalized:
params.g.normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
params.g.normalizationRegion = 'zscore'; % 'none', 'zscore'

%----Compute enrichment results
params.c.structFilter = 'all';
results.mouse_all = geneEnrichmentDistance(params);

% Let's try to visualize gene categories with different types of scores:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);

categoryIndex = 9;

isThisCat = ismember(geneInfo.entrez_id,results.annotations{categoryIndex});
geneDataSub = geneData(:,isThisCat);

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
title(results.GOName{categoryIndex})

% Spatial projection:
subplot(1,3,2); axis('equal')
% nanmean(geneDataSub,2);
gMeanFlipped = flippedMean(geneDataSub);
gMean = BF_NormalizeMatrix(gMeanFlipped,'mixedSigmoid');
distMat = GiveMeDistanceMatrix(params.humanOrMouse);
coOrdsXY = mdscale(distMat,2);
scatter(coOrdsXY(:,1),coOrdsXY(:,2),25,gMean,'filled')
colormap([flipud(BF_getcmap('blues',9));1,1,1;BF_getcmap('reds',9)])
xlabel('spatialAxis_1')
ylabel('spatialAxis_2')
colorbar

% Histogram of individual gene correlations:
ax = subplot(1,3,3); hold on
isUpperDiag = triu(true(size(distMat)),+1);
individualRhos = zeros(size(geneDataSub,2),1);
for j = 1:size(geneDataSub,2)
    GCCj = geneDataSub(:,j)*geneDataSub(:,j)';
    individualRhos(j) = corr(distMat(isUpperDiag),GCCj(isUpperDiag));
end
histogram(individualRhos)
GCCmean = gMean*gMean';
rhoMean = corr(distMat(isUpperDiag),GCCmean(isUpperDiag))
plot(ones(2,1)*rhoMean,ax.YLim,'r')
plot(-ones(2,1)*results.meanScore(categoryIndex),ax.YLim,'g')
xlabel('\rho (GCC-distance)')
% plot(distMat(isUpperDiag),GCC(isUpperDiag),'.k')
% ylabel('Mean expression GCC')
% title(sprintf('rho = %.3f',rho))

%-------------------------------------------------------------------------------
params.c.structFilter = 'isocortex';
results.mouse_ctx = geneEnrichmentDistance(params);

%-------------------------------------------------------------------------------
% Human:
params = GiveMeDefaultParams('human');
params.g.whatParcellation = 'HCP';
params.g.normalizationInternal = 'robustSigmoid';
params.g.normalizationGene = 'mixedSigmoid'; % 'none', 'mixedSigmoid'
params.g.normalizationRegion = 'none'; % 'none', 'zscore'

results.human_HCP = geneEnrichmentDistance(params);

%-------------------------------------------------------------------------------
% Visualize together:
thresholdSig = 0.05;
PlotEnrichmentTables(results,thresholdSig);
title('Enrichment by expression variance across the brain')
