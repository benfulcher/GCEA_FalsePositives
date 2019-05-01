%-------------------------------------------------------------------------------
% MouseRandom
%-------------------------------------------------------------------------------

% Reproduce a massive calculation
computeMode = false;

%===============================================================================
% Within-category correlation
%===============================================================================
if computeMode
    numSamples = 20000;
    params = GiveMeDefaultParams('mouse');
    mouseIntra = IntraCorrelationByCategory(params,'geneShuffle',numSamples,'VE1',true);
end

%===============================================================================
% Import intra-data and random data
%===============================================================================
results = struct();
resultsIntra = load('mouseIntra_geneShuffle_20k.mat');
results.intra = resultsIntra.mouseIntra;
results.random = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','');
results.randomNull = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','independentSpatialShuffle');

[rowVectorResults,GOTerms,allGOIDs,tableNames] = CombineTables(results,'mouse',{'pValZCorr','sumUnderSig','sumUnderSig'});

%-------------------------------------------------------------------------------
f = figure('color','w');
plot(rowVectorResults(1,:),rowVectorResults(2,:)/1e4,'.k')
xlabel(sprintf('%s-meanIntraClassCorrelation',tableNames{1}))
ylabel(sprintf('%s-propSignificantRandomPermutation',tableNames{2}))
axis('square')
hold on
plot(rowVectorResults(1,:),rowVectorResults(3,:)/1e4,'.r')

f = figure('color','w');
hold('on');
BF_PlotQuantiles(rowVectorResults(1,:),rowVectorResults(2,:)/1e4,20,false,false,'k');
BF_PlotQuantiles(rowVectorResults(1,:),rowVectorResults(3,:)/1e4,20,false,false,'r');
ylabel('Probability of significance under random null')
xlabel('P-value of intra-category correlation')

%-------------------------------------------------------------------------------
results = struct();
results.randomReal = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','');
results.randomSpatialCoord = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','coordinatedSpatialShuffle');
results.randomSpatialInd = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','independentSpatialShuffle');
