function MouseRandom(whatSpecies,whatShuffle)
%-------------------------------------------------------------------------------
% MouseRandom
%-------------------------------------------------------------------------------
if nargin < 1
    whatSpecies = 'mouse';
end
if nargin < 2
    whatShuffle = 'geneShuffle'; % 'geneShuffle', 'independentSpatialShuffle'
end

%-------------------------------------------------------------------------------
% Reproduce a massive calculation
computeMode = false;
numNullSamples_VE1 = 20000; % (Intra_*_VE1_20000.mat)
numNullSamples_surrogate = 10000; % (SurrogateGOTables_10000_*.mat)

%===============================================================================
% Within-category correlation
%===============================================================================
if computeMode
    params = GiveMeDefaultParams('mouse');
    mouseIntra = IntraCorrelationByCategory(params,whatShuffle,numNullSamples_VE1,'VE1',true);
end

%===============================================================================
% Import intra-data and surrogate random data
%===============================================================================
results = struct();
resultsIntra = load(sprintf('Intra_%s_%s_VE1_%u.mat',whatSpecies,whatShuffle,numNullSamples_VE1));
results.intra = resultsIntra.resultsTable;
results.randomReal = SurrogateEnrichmentProcess('mouse',numNullSamples_surrogate,'randomUniform','');
% results.randomSpatialCoord = SurrogateEnrichmentProcess('mouse',numNullSamples_surrogate,'randomUniform','coordinatedSpatialShuffle');
results.randomSpatialInd = SurrogateEnrichmentProcess('mouse',numNullSamples_surrogate,'randomUniform','independentSpatialShuffle');

[rowVectorResults,GOTerms,allGOIDs,tableNames] = CombineTables(results,'mouse',{'pValZCorr','sumUnderSig','sumUnderSig'});

%-------------------------------------------------------------------------------
f = figure('color','w');
hold('on')
pRandomNull = plot(rowVectorResults(1,:),rowVectorResults(2,:)/numNullSamples_surrogate,'.k');
pIntraNull = plot(rowVectorResults(1,:),rowVectorResults(3,:)/numNullSamples_surrogate,'.r');
xlabel(sprintf('%s-IntraClassVE1',tableNames{1}));
ylabel(sprintf('%s-propSignificantRandomPermutation',tableNames{2}));
axis('square')
legend('intraCorr','randomChance')

%-------------------------------------------------------------------------------
f = figure('color','w');
hold('on');
numBins = 10;
BF_PlotQuantiles(rowVectorResults(1,:),rowVectorResults(2,:)/numNullSamples_surrogate,numBins,false,false,'k',false);
BF_PlotQuantiles(rowVectorResults(1,:),rowVectorResults(3,:)/numNullSamples_surrogate,numBins,false,false,[27,158,119]/255,false);
ylabel('Probability of p < 0.05 (random null)')
xlabel('p-value of intra-category VE1')
f.Position = [1000        1117         326         221];


end
