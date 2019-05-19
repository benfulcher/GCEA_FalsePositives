
results = struct();

% Interspecies correspondence of VE1 p-values:
resultsIntraMouse = load('Intra_mouse_geneShuffle_VE1_20000.mat');
results.intraMouse = resultsIntraMouse.resultsTable;

resultsIntraHuman = load('Intra_human_geneShuffle_VE1_20000.mat');
results.intraHuman = resultsIntraHuman.resultsTable;

[rowVectorResults,GOTerms,allGOIDs,tableNames] = CombineTables(results,'mouse',{'intracorr_VE1','intracorr_VE1'});


f = figure('color','w');
plot(log10(rowVectorResults(1,:)),log10(rowVectorResults(2,:)),'.k')
xlabel(tableNames{1})
ylabel(tableNames{2})

[r,p] = corr(rowVectorResults(1,:)',rowVectorResults(2,:)','type','Spearman',...
                                    'rows','pairwise')
