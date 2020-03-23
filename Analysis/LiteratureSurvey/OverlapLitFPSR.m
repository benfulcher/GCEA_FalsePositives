function OverlapLitFPSR(whatSpecies,makeFigure)
% Understand the overlap between high FPSR categories and those reported in the
% literature
%-------------------------------------------------------------------------------
if nargin < 1
    whatSpecies = 'human';
end
if nargin < 2
    makeFigure = false;
end
params = GiveMeDefaultParams(whatSpecies);

%-------------------------------------------------------------------------------
% Retrieve information about how literature results are distributed across categories:
LitTable = MakeLiteratureTable(whatSpecies,params.e.sigThresh);

%-------------------------------------------------------------------------------
% Now we'll get the FPSR data:
GOTable_FPSR = SurrogateEnrichmentProcess(params,false);

%-------------------------------------------------------------------------------
% LITERATURE categories: LitTable
% If at least X studies have quoted it at pCorr < 0.05, then we label it as a "reported category":
GOID_reported_1 = LitTable.GOID(LitTable.numStudies==1);
GOID_reported_2 = LitTable.GOID(LitTable.numStudies==2);
GOID_reported_3plus = LitTable.GOID(LitTable.numStudies>=3);
GOID_reported_2plus = LitTable.GOID(LitTable.numStudies>=2);

% FPSR categories, GOTable_FPSR
% Split into reported and unreported categories:
% Was a category reported anywhere in the literature?
GOTable_notReported = GOTable_FPSR(~ismember(GOTable_FPSR.GOID,LitTable.GOID),:);
GOTable_reported1 = GOTable_FPSR(ismember(GOTable_FPSR.GOID,GOID_reported_1),:);
GOTable_reported2 = GOTable_FPSR(ismember(GOTable_FPSR.GOID,GOID_reported_2),:);
GOTable_reported_3plus = GOTable_FPSR(ismember(GOTable_FPSR.GOID,GOID_reported_3plus),:);
GOTable_reported_2plus = GOTable_FPSR(ismember(GOTable_FPSR.GOID,GOID_reported_2plus),:);
GOTable_notreported_2plus = GOTable_FPSR(~ismember(GOTable_FPSR.GOID,GOID_reported_2plus),:);

%-------------------------------------------------------------------------------
% PLOT: Distributions

% FPSR_reported = GOTable_reported_2plus.sumUnderSig;
% FPSR_notReported = GOTable_notreported_2plus.sumUnderSig;
% BF_JitteredParallelScatter({FPSR_reported,FPSR_notReported},true,true,true);

%-------------------------------------------------------------------------------
% PLOT: Proportion by bins

switch whatSpecies
case 'mouse'
    numBins = 5;
    theFields = [2,5]; % 1, or 2+
    legendEntries = {'1 study','2+ studies'};
case 'human'
    numBins = 5;
    theFields = [3,4]; % 2, or 3+
    legendEntries = {'2 studies','3+ studies'};
end
binEdges = linspace(0,max(GOTable_FPSR.sumUnderSig)+eps,numBins+1);
FPSR_repOrNot = zeros(numBins,5);
isInBin = @(x,Gtable) sum(Gtable.sumUnderSig >= binEdges(x) & Gtable.sumUnderSig < binEdges(x+1));
for i = 1:numBins
    FPSR_repOrNot(i,1) = isInBin(i,GOTable_notReported);
    FPSR_repOrNot(i,2) = isInBin(i,GOTable_reported1);
    FPSR_repOrNot(i,3) = isInBin(i,GOTable_reported2);
    FPSR_repOrNot(i,4) = isInBin(i,GOTable_reported_3plus);
    FPSR_repOrNot(i,5) = isInBin(i,GOTable_reported_2plus);
end
propReportedNorm = FPSR_repOrNot./sum(FPSR_repOrNot,2);
binMeans = mean([binEdges(1:end-1);binEdges(2:end)]);
if makeFigure
    f = figure('color','w');
end
bh = bar(100*binMeans/10000,100*propReportedNorm(:,theFields),'stacked','BarWidth',0.95);
xlabel('CFPR (%)')
ylabel('GO categories reported as significant (%)')
legend(legendEntries,'Location','northwest')
title(whatSpecies)
colors = [97,184,255; 255,163,0]/255;
bh(1).FaceColor = colors(1,:);
bh(2).FaceColor = colors(2,:);
ax = gca;
ax.XLim = [0,30];
% f.Position = [1058         972         460         219];
f.Position = [1000        1121         273         136];
% bar(binMeans/10000,FPSR_repOrNot,'stacked')

% Save:
if makeFigure
    fileName = fullfile('OutputPlots',sprintf('CFPR_literature_histogram_%s.svg',whatSpecies));
    saveas(f,fileName,'svg')
    fprintf(1,'Saved to %s\n',fileName);
end

end
