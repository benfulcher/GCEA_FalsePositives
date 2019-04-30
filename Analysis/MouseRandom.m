%-------------------------------------------------------------------------------
% MouseRandom
%-------------------------------------------------------------------------------

%===============================================================================
% Within-category correlation
%===============================================================================
numSamples = 100000;
params = GiveMeDefaultParams('mouse');
mouseIntra = IntraCorrelationByCategory(params,'geneShuffle',numSamples);
fileOut = fullfile('DataOutputs','mouseIntra_geneShuffle_20k.mat')
save(fileOut,'mouseIntra','params','numSamples');

%===============================================================================
% Import intra-data and random data
%===============================================================================
results = struct();
resultsIntra = load(fullfile('DataOutputs','mouseIntra_geneShuffle_20k.mat'));
results.intra = resultsIntra.mouseIntra;
results.random = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','');

[rowVectorResults,GOTerms,allGOIDs,tableNames] = CombineTables(results,'mouse',{'pValZ','pValCorr'});

%-------------------------------------------------------------------------------
f = figure('color','w');
plot(rowVectorResults(1,:),rowVectorResults(2,:),'.k')
xlabel(tableNames{1})
ylabel(tableNames{2})

%-------------------------------------------------------------------------------
results = struct();
results.randomReal = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','');
results.randomSpatialCoord = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','coordinatedSpatialShuffle');
results.randomSpatialInd = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','independentSpatialShuffle');
