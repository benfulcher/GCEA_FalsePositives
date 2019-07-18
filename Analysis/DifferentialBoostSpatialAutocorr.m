% Does adding spatial autocorrelation to phenotypes boost the correlations to
% GO categories containing genes with high spatial autocorrelation?

numNullSamples_surrogate = 10000;
whatSpecies = 'mouse'; structFilter = 'all';
% whatSpecies = 'human'; structFilter = 'cortex';
params = GiveMeDefaultParams(whatSpecies,structFilter);


% Get spatial autocorrelation scores per GO category:
params.g.normalizationGene = 'zscore';
params.g.normalizationRegion = 'zscore';
results = geneEnrichmentDistance(params);

% Load FPSR (random)
GOTable_FPSR_mouseRandom = SurrogateEnrichmentProcess(whatSpecies,numNullSamples_surrogate,'randomUniform','');
% Load FPSR (spatial autocorrelation)
GOTable_FPSR_mouseAC = SurrogateEnrichmentProcess(whatSpecies,numNullSamples_surrogate,'spatialLag','');

% Differential FPSR:
% Let's combine/annotate:
[~,ia,ib] = intersect(GOTable_FPSR_mouseRandom.GOID,GOTable_FPSR_mouseAC.GOID);
GOTable_FPSR_mouseRandom = GOTable_FPSR_mouseRandom(ia,:);
GOTable_FPSR_mouseAC = GOTable_FPSR_mouseAC(ib,:);

GOTable_FPSR_mouseRandom.sumUnderSig_AC = GOTable_FPSR_mouseAC.sumUnderSig;
GOTable_FPSR_mouseRandom.absDiffFPSR = (GOTable_FPSR_mouseRandom.sumUnderSig_AC - GOTable_FPSR_mouseRandom.sumUnderSig)/100;
GOTable_FPSR_mouseRandom.relDiffFPSR = 100*(GOTable_FPSR_mouseRandom.sumUnderSig_AC - GOTable_FPSR_mouseRandom.sumUnderSig)./GOTable_FPSR_mouseRandom.sumUnderSig;

%-------------------------------------------------------------------------------
f = figure('color','w');
histogram(GOTable_FPSR_mouseRandom.absDiffFPSR)
xlabel('Absolute FPSR (spatialLag - random) (%)')
ylabel('Frequency')
title(whatSpecies)

%-------------------------------------------------------------------------------
% Does the level of spatial autocorrelation of genes in a category capture the change?:
[~,ia,ib] = intersect(GOTable_FPSR_mouseRandom.GOID,results.GOID);
GOTableCombined = GOTable_FPSR_mouseRandom(ia,:);
% annotate the meanScore
GOTableCombined.meanACscore = results.meanScore(ib);
GOTableCombined.pValZCorr = results.pValZCorr(ib);

%-------------------------------------------------------------------------------
f = figure('color','w');
plot(GOTableCombined.meanACscore,GOTableCombined.relDiffFPSR,'.k')
xlabel('mean spatial AC score')
ylabel('Relative FPSR (spatialLag - random) (%)')

%-------------------------------------------------------------------------------
f = figure('color','w');
numBins = 20;
theColor = BF_getcmap('spectral',3,1);
BF_PlotQuantiles(GOTableCombined.meanACscore,GOTableCombined.relDiffFPSR,numBins,false,false,'k',false);
xlabel('Mean autocorrelation score')
ylabel('Relative change FPSR(AC) - FPSR(random) (%)')
title(whatSpecies)


%===============================================================================
%===============================================================================
% Intra-category correlation:
numNullSamples_intraCorr = 20000; % (Intra_*_*_20000.mat)
numNullSamples_surrogate = 10000; % (SurrogateGOTables_10000_*.mat)
whatIntraStat = 'raw';
whatShuffle = 'geneShuffle';
fileNameIn = sprintf('Intra_%s_%s_%s_%u.mat',whatSpecies,whatShuffle,whatIntraStat,numNullSamples_intraCorr);
resultsIntra = load(fileNameIn);
resultsIntra = resultsIntra.resultsTable;



[~,ia,ib] = intersect(GOTableCombined.GOID,resultsIntra.GOID);
GOTableCombinedAgain = GOTableCombined(ia,:);
GOTableCombinedAgain.intracorr_raw = resultsIntra.intracorr_raw(ib);

f = figure('color','w'); hold on
xData = GOTableCombinedAgain.intracorr_raw;
yData = GOTableCombinedAgain.meanACscore;
colorData = GOTableCombinedAgain.sumUnderSig;
backGround = (colorData < 500);
plot(xData(backGround),yData(backGround),'.k');
scatter(xData(~backGround),yData(~backGround),...
                        50,colorData(~backGround),'filled','MarkerEdgeColor','k')
title(whatSpecies)
xlabel('intra-category coexpression')
ylabel('mean spatial correlation')
colorbar
colormap(flipud(BF_getcmap('redyellowblue',8)));
