

whatSpecies = 'mouse';

nullGOTables = struct();
nullGOTables.spatialRandom = SurrogateEnrichmentProcess(whatSpecies,10000,'randomUniform','');
% nullGOTables.coordSpatialRandom = SurrogateEnrichmentProcess(whatSpecies,10000,'randomUniform','coordinatedSpatialShuffle');
nullGOTables.indSpatialRandom = SurrogateEnrichmentProcess(whatSpecies,10000,'randomUniform','independentSpatialShuffle');
nullGOTables.spatialLag = SurrogateEnrichmentProcess(whatSpecies,5000,'spatialLag','');

%-------------------------------------------------------------------------------
% Extract sum under significant values:
sumUnderSigValues = structfun(@(x)x.sumUnderSig/10000,nullGOTables,'UniformOutput',false);
sumUnderSigValues.spatialLag = sumUnderSigValues.spatialLag*2;% (half as many null samples)

%-------------------------------------------------------------------------------
% Distribution of sumUnderSig:
fields = fieldnames(nullGOTables);
sumUnderSigCell = struct2cell(sumUnderSigValues);
sumUnderSigCell = sumUnderSigCell([3,1,2])
% f = figure('color','w'); hold('on')
% histogram(sumUnderSigCell{3},'normalization','probability')
% histogram(sumUnderSigCell{1},'normalization','probability')
% histogram(sumUnderSigCell{2},'normalization','probability')
BF_JitteredParallelScatter(sumUnderSigCell)
ax = gca();
ax.XTick = 1:4;
ax.XTickLabel = fields([3,1,2]);
ylabel('Proportion significant under random null')
title('Distribution over GO categories')
f = gcf();
f.Position = [1000        1121         273         217];
