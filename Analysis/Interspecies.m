
% Parameters:
whatShuffle = 'independentSpatialShuffle'; % 'geneShuffle'
numSamples = 20000;
theIntraCategoryStat = 'VE1';
%-------------------------------------------------------------------------------

results = struct();

% Interspecies correspondence of VE1 p-values:
resultsIntraMouse = load(sprintf('Intra_mouse_%s_%s_%u.mat',whatShuffle,theIntraCategoryStat,numSamples));
results.intraMouse = resultsIntraMouse.resultsTable;

resultsIntraHuman = load(sprintf('Intra_human_%s_%s_%u.mat',whatShuffle,theIntraCategoryStat,numSamples));
results.intraHuman = resultsIntraHuman.resultsTable;

% [rowVectorResults,allGOIDs,tableNames] = CombineTables(results,'mouse',{'pValZ','pValZ'});
[rowVectorResults,allGOIDs,tableNames] = CombineTables(results,'mouse',{'intracorr_VE1','intracorr_VE1'});


f = figure('color','w');
plot(rowVectorResults(1,:),rowVectorResults(2,:),'.k')
% plot(log10(rowVectorResults(1,:)),log10(rowVectorResults(2,:)),'.k')
xlabel(sprintf('%s-%s',tableNames{1},theIntraCategoryStat))
ylabel(sprintf('%s-%s',tableNames{2},theIntraCategoryStat))

[rho,p] = corr(rowVectorResults(1,:)',rowVectorResults(2,:)','type','Spearman',...
                                    'rows','pairwise')
