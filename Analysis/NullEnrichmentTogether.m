

whatSpecies = 'mouse';

nullGOTables = struct();
nullGOTables.spatialRandom = SurrogateEnrichmentProcess(whatSpecies,10000,'randomUniform','');
nullGOTables.coordSpatialRandom = SurrogateEnrichmentProcess(whatSpecies,10000,'randomUniform','coordinatedSpatialShuffle');
nullGOTables.indSpatialRandom = SurrogateEnrichmentProcess(whatSpecies,10000,'randomUniform','independentSpatialShuffle');
nullGOTables.spatialLag = SurrogateEnrichmentProcess(whatSpecies,5000,'spatialLag','');

%-------------------------------------------------------------------------------
% Extract sum under significant values:
sumUnderSigValues = structfun(@(x)x.sumUnderSig/10000,nullGOTables,'UniformOutput',false);
sumUnderSigValues.spatialLag = sumUnderSigValues.spatialLag*2;% (half as many null samples)

%-------------------------------------------------------------------------------
% Distribution of sumUnderSig:
fields = fieldnames(nullGOTables);
BF_JitteredParallelScatter(struct2cell(sumUnderSigValues))
ax = gca();
ax.XTick = 1:4;
ax.XTickLabel = fields;
ylabel('Proportion significant under random null')
title('Distribution over GO categories')
