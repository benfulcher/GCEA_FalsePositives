function MouseRandom(whatSpecies,whatShuffle,whatIntraStat)
%-------------------------------------------------------------------------------
% MouseRandom
%-------------------------------------------------------------------------------
if nargin < 1
    whatSpecies = 'human';
end
if nargin < 2
    whatShuffle = 'geneShuffle'; % 'geneShuffle', 'independentSpatialShuffle'
end
if nargin < 3
    whatIntraStat = 'raw';
end

%-------------------------------------------------------------------------------
computeMode = false; % Reproduce a massive calculation
numNullSamples_intraCorr = 20000; % (Intra_*_*_20000.mat)
numNullSamples_surrogate = 10000; % (SurrogateGOTables_10000_*.mat)

%===============================================================================
% Import (or compute) intra-category correlation data and surrogate random data
%===============================================================================
results = struct();
params = GiveMeDefaultParams(whatSpecies);
if computeMode
    resultsIntra = IntraCorrelationByCategory(params,whatShuffle,numNullSamples_intraCorr,whatIntraStat,true);
else
    fileNameIn = sprintf('Intra_%s_%s_%s_%u.mat',whatSpecies,whatShuffle,whatIntraStat,numNullSamples_intraCorr);
    resultsIntra = load(fileNameIn);
    fprintf(1,'Importing intra-category coexpression data from %s\n',fileNameIn);
end
results.intra = resultsIntra.resultsTable;
results.randomReal = SurrogateEnrichmentProcess(whatSpecies,numNullSamples_surrogate,'randomUniform','');
% results.randomSpatialInd = SurrogateEnrichmentProcess(whatSpecies,numNullSamples_surrogate,'randomUniform','independentSpatialShuffle');
% [rowVectorResults,allGOIDs,tableNames] = CombineTables(results,whatSpecies,{'intracorr_raw','sumUnderSig','sumUnderSig'});

[~,ia,ib] = intersect(results.intra.GOID,results.randomReal.GOID);
GOTableCombined = results.intra(ia,:);
GOTableCombined.FPSR_random = results.randomReal.sumUnderSig(ib)/numNullSamples_surrogate;

%-------------------------------------------------------------------------------
% What about as a function of category size:
numBins = 1;
binEdges = linspace(min(GOTableCombined.size),max(GOTableCombined.size)+eps,numBins+1);
isInBin = @(x) (GOTableCombined.size>=binEdges(x) & GOTableCombined.size < binEdges(x+1));
binIndex = arrayfun(@(x)isInBin(x),1:numBins,'UniformOutput',false);

%-------------------------------------------------------------------------------
f = figure('color','w');
theXfield = 'intracorr_raw';
theYfield = 'FPSR_random';
hold('on')
if numBins==1
    plot(GOTableCombined.(theXfield),GOTableCombined.(theYfield),'.k');
else
    colors = BF_getcmap('spectral',numBins,1);
    ph = cell(numBins,1);
    for k = 1:numBins
        ph{k} = plot(GOTableCombined.(theXfield)(binIndex{k}),GOTableCombined.(theYfield)(binIndex{k}),'.','color',colors{k});
    end
    legend([ph{:}],arrayfun(@(x)sprintf('Bin%u',x),1:numBins,'UniformOutput',false))
end
xlabel(sprintf('%s-IntraClass%s',tableNames{1},whatIntraStat));
ylabel(sprintf('%s-propSignificantRandomPermutation',tableNames{2}));
axis('square')

%-------------------------------------------------------------------------------
f = figure('color','w'); hold('on');
numQuantileBins = 10;
theColor = BF_getcmap('spectral',3,1);
if numBins==1
    BF_PlotQuantiles(GOTableCombined.(theXfield),GOTableCombined.(theYfield),...
            numQuantileBins,false,false,theColor{1},false);
else
    for k = 1:numBins
        BF_PlotQuantiles(GOTableCombined.(theXfield)(binIndex{k}),GOTableCombined.(theYfield)(binIndex{k}),...
                numQuantileBins,false,false,colors{k},false);
    end
end
% BF_PlotQuantiles(rowVectorResults(1,:),rowVectorResults(3,:)/numNullSamples_surrogate,numBins,false,false,[27,158,119]/255,false);
ylabel('FPSR (random null)')
xlabel(sprintf('intra-category %s',whatIntraStat))
f.Position = [1000        1117         326         221];
title(whatSpecies)

%-------------------------------------------------------------------------------
% Is intra-category correlation related to that category's null width?
% Estimate null width first using FPSR

% Take category sizes from results.randomReal:
[rowVectorResults,allGOIDs,tableNames] = CombineTables(results,whatSpecies,{'intracorr_VE1','sumUnderSig','sumUnderSig'});
[allGOIDsMatched,ia,ib] = intersect(allGOIDs,results.randomReal.GOID,'stable');
categorySizesMatched = results.randomReal.size(ib);
intracorr_VE1 = rowVectorResults(1,ia)';
sumUnderSig = rowVectorResults(2,ia)';

% Match to estimates of null width under random spatial maps:
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

% Compute correlations with VE1 in binned category sizes:
numBins = 10;
numThresholds = numBins + 1;
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

fprintf('Average over %u bins of category size: Spearman correlation of VE1 with FPSR: %.2g\n',...
                        numBins,mean(corrInSizeFPSR(:,1)));
fprintf('Average over %u bins of category size: Spearman correlation of VE1 with random-map null width: %.2g\n',...
                        numBins,nanmean(corrInSizeNullWidth(:,1)));

%-------------------------------------------------------------------------------
f = figure('color','w'); hold('on')
xlabel('Category size bins')
ylabel('Mean FPSR in bin')
title(whatSpecies)
for k = 1:numBins
    plot(xThresholds(k:k+1),ones(2,1)*corrInSizeFPSR(k),'LineStyle','-','LineWidth',2,'Color','k')
end


%-------------------------------------------------------------------------------
% What about the distinctiveness of the genes in the category
% GOTable = AnnotateDistinctivenessScore(params);
% [rowVectorResults,allGOIDs,tableNames] = CombineTables(results,whatSpecies,{'intracorr_VE1','sumUnderSig','sumUnderSig'});



end
