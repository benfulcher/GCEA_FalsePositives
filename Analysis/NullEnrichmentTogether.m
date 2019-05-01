


whatSpecies = 'mouse';

nullGOTables = struct();
nullGOTables.spatialRandom = SurrogateEnrichmentProcess(whatSpecies,10000,'randomUniform','');
nullGOTables.coordSpatialRandom = SurrogateEnrichmentProcess(whatSpecies,10000,'randomUniform','coordinatedSpatialShuffle');
nullGOTables.indSpatialRandom = SurrogateEnrichmentProcess(whatSpecies,10000,'randomUniform','independentSpatialShuffle');
% nullGOTables.spatialLag = SurrogateEnrichmentProcess(whatSpecies,1000,'spatialLag','');

%-------------------------------------------------------------------------------
% Distribution of sumUnderSig:
sumUnderSigValues = structfun(@(x)x.sumUnderSig/10000,nullGOTables,'UniformOutput',false);
fields = fieldnames(nullGOTables);

BF_JitteredParallelScatter(struct2cell(sumUnderSigValues))
ax = gca();
ax.XTick = 1:3;
ax.XTickLabel = fields;
ylabel('Proportion significant under random null')
title('Distribution over GO categories')
