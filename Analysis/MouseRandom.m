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
whatShuffle = 'independentSpatialShuffle';
%-------------------------------------------------------------------------------
% Reproduce a massive calculation
computeMode = false;
numNullSamples_VE1 = 20000; % (Intra_*_VE1_20000.mat)
numNullSamples_surrogate = 10000; % (SurrogateGOTables_10000_*.mat)

%===============================================================================
% Import (or compute) intra-category correlation data and surrogate random data
%===============================================================================
results = struct();
if computeMode
    params = GiveMeDefaultParams(whatSpecies);
    resultsIntra = IntraCorrelationByCategory(params,whatShuffle,numNullSamples_VE1,'VE1',true);
else
    resultsIntra = load(sprintf('Intra_%s_%s_VE1_%u.mat',whatSpecies,whatShuffle,numNullSamples_VE1));
end
results.intra = resultsIntra.resultsTable;
results.randomReal = SurrogateEnrichmentProcess(whatSpecies,numNullSamples_surrogate,'randomUniform','');
results.randomSpatialInd = SurrogateEnrichmentProcess(whatSpecies,numNullSamples_surrogate,'randomUniform','independentSpatialShuffle');

[rowVectorResults,allGOIDs,tableNames] = CombineTables(results,whatSpecies,{'pValZCorr','sumUnderSig','sumUnderSig'});

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

%-------------------------------------------------------------------------------
% Is intra-category correlation related to that category's null width?
% Estimate null width first using FPSR

% Take category sizes from results.randomReal:
[rowVectorResults,allGOIDs,tableNames] = CombineTables(results,whatSpecies,{'intracorr_VE1','sumUnderSig','sumUnderSig'});
[allGOIDsMatched,ia,ib] = intersect(allGOIDs,results.randomReal.GOID,'stable');
categorySizesMatched = results.randomReal.size(ib);
intracorr_VE1 = rowVectorResults(1,ia)';
sumUnderSig = rowVectorResults(2,ia)';

% Match to better estimates (null width):
numNullSamples = 20000;
whatSurrogate = 'randomMap';
whatCorr = 'Spearman';
theDataFile = sprintf('RandomNull_%u_%s-%s_%s_%s_mean.mat',numNullSamples,whatSpecies,...
                            resultsIntra.params.g.structFilter,whatSurrogate,whatCorr);
fprintf(1,'Loading in precomputed null data from ''%s''\n',theDataFile);
load(theDataFile,'GOTable');
randomMapNullWidth = nan(length(allGOIDsMatched),1);
for i = 1:length(allGOIDsMatched)
    whatCategory = find(GOTable.GOID==allGOIDsMatched(i));
    if ~isempty(whatCategory)
        randomMapNullWidth(i) = var(GOTable.categoryScores{whatCategory});
    end
end

% Bin the category sizes:
numThresholds = 11;
numBins = numThresholds-1;
% xThresholds = arrayfun(@(x)quantile(categorySizesMatched,x),linspace(0,1,numThresholds));
xThresholds = round(linspace(min(categorySizesMatched),max(categorySizesMatched),numThresholds));
xThresholds(end) = xThresholds(end) + eps; % make sure all data is included in final bin
corrInSizeFPSR = zeros(numBins,2);
corrInSizeNullWidth = zeros(numBins,2);
for i = 1:numBins
    isInThreshold = (categorySizesMatched>=xThresholds(i) & categorySizesMatched < xThresholds(i+1));
    [corrInSizeFPSR(i,1),corrInSizeFPSR(i,2)] = corr(intracorr_VE1(isInThreshold),sumUnderSig(isInThreshold),'type','Spearman','rows','pairwise');
    [corrInSizeNullWidth(i,1),corrInSizeNullWidth(i,2)] = corr(randomMapNullWidth(isInThreshold),sumUnderSig(isInThreshold),'type','Spearman','rows','pairwise');
end
[corrIgnoreSizeFPSR,p] = corr(intracorr_VE1,sumUnderSig,'type','Spearman','rows','pairwise');
[corrIgnoreSizeNullWidth,p] = corr(randomMapNullWidth,sumUnderSig,'type','Spearman','rows','pairwise');

f = figure('color','w'); hold('on')
for k = 1:numBins
    plot(xThresholds(k:k+1),ones(2,1)*corrInSizeFPSR(k),'LineStyle','-','LineWidth',2,'Color','k')
end


end
