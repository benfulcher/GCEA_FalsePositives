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
resultsIntra = load('Intra_mouse_geneShuffle_VE1_20000.mat');
results.intra = resultsIntra.resultsTable;
results.random = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','');
results.randomNull = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','independentSpatialShuffle');

[rowVectorResults,GOTerms,allGOIDs,tableNames] = CombineTables(results,'mouse',{'pValZCorr','sumUnderSig','sumUnderSig'});

%-------------------------------------------------------------------------------
f = figure('color','w');
hold('on')
pRandomNull = plot(rowVectorResults(1,:),rowVectorResults(2,:)/1e4,'.k');
pIntraNull = plot(rowVectorResults(1,:),rowVectorResults(3,:)/1e4,'.r');
xlabel(sprintf('%s-IntraClassVE1',tableNames{1}));
ylabel(sprintf('%s-propSignificantRandomPermutation',tableNames{2}));
axis('square')
legend('intraCorr','randomChance')

%-------------------------------------------------------------------------------
f = figure('color','w');
hold('on');
numBins = 10;
BF_PlotQuantiles(rowVectorResults(1,:),rowVectorResults(2,:)/1e4,numBins,false,false,'k',false);
BF_PlotQuantiles(rowVectorResults(1,:),rowVectorResults(3,:)/1e4,numBins,false,false,[27,158,119]/255,false);
ylabel('Probability of p < 0.05 (random null)')
xlabel('p-value of intra-category VE1')
f.Position = [1000        1117         326         221];

%-------------------------------------------------------------------------------
results = struct();
results.randomReal = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','');
results.randomSpatialCoord = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','coordinatedSpatialShuffle');
results.randomSpatialInd = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','independentSpatialShuffle');
