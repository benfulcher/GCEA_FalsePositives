% Set parameters:
whatSpecies = 'human';
pValCorrThreshold = 0.05;

%-------------------------------------------------------------------------------
% Retrieve information about how literature results are distributed across categories:
LitTable = MakeLiteratureTable(whatSpecies,pValCorrThreshold);

%-------------------------------------------------------------------------------
% Now we'll get the FPSR data:
numNullSamples_surrogate = 10000;
GOTable_FPSR = SurrogateEnrichmentProcess(whatSpecies,numNullSamples_surrogate,'randomUniform','');

% If at least X studies have quoted it at pCorr < 0.05, then we label it as a "reported category":
GOID_reported_1 = LitTable.GOID(LitTable.numStudies==1);
GOID_reported_2 = LitTable.GOID(LitTable.numStudies==2);
GOID_reported_3plus = LitTable.GOID(LitTable.numStudies>=3);
GOID_reported_2plus = LitTable.GOID(LitTable.numStudies>=2);

% Split into reported and unreported categories:
GOTable_notReported = GOTable_FPSR(~ismember(GOTable_FPSR.GOID,LitTable.GOID),:);
GOTable_reported1 = GOTable_FPSR(ismember(GOTable_FPSR.GOID,GOID_reported_1),:);
GOTable_reported2 = GOTable_FPSR(ismember(GOTable_FPSR.GOID,GOID_reported_2),:);
GOTable_reported_3plus = GOTable_FPSR(ismember(GOTable_FPSR.GOID,GOID_reported_3plus),:);

GOTable_reported_2plus = GOTable_FPSR(ismember(GOTable_FPSR.GOID,GOID_reported_2plus),:);
GOTable_notreported_2plus = GOTable_FPSR(~ismember(GOTable_FPSR.GOID,GOID_reported_2plus),:);

% Distributions:
FPSR_reported = GOTable_reported_2plus.sumUnderSig;
FPSR_notReported = GOTable_notreported_2plus.sumUnderSig;
BF_JitteredParallelScatter({FPSR_reported,FPSR_notReported},true,true,true);

% Propotion by bins:
switch whatSpecies
case 'mouse'
    numBins = 5;
    theFields = [2,5]; % 1, or 2+
    legendEntries = {'1 study','2+ studies'};
case 'human'
    numBins = 5;
    theFields = [3,4]; % 1, or 2+
    legendEntries = {'2 studies','3+ studies'};
end
binEdges = linspace(0,max(GOTable_FPSR.sumUnderSig)+eps,numBins+1);
FPSR_repOrNot = zeros(numBins,5);
isInBin = @(x,Gtable) sum(Gtable.sumUnderSig>=binEdges(x) & Gtable.sumUnderSig < binEdges(x+1));
for i = 1:numBins
    FPSR_repOrNot(i,1) = isInBin(i,GOTable_notReported);
    FPSR_repOrNot(i,2) = isInBin(i,GOTable_reported1);
    FPSR_repOrNot(i,3) = isInBin(i,GOTable_reported2);
    FPSR_repOrNot(i,4) = isInBin(i,GOTable_reported_3plus);
    FPSR_repOrNot(i,5) = isInBin(i,GOTable_reported_2plus);
end
propReportedNorm = FPSR_repOrNot./sum(FPSR_repOrNot,2);
binMeans = mean([binEdges(1:end-1);binEdges(2:end)]);
f = figure('color','w');
bh = bar(100*binMeans/10000,100*propReportedNorm(:,theFields),'stacked','BarWidth',0.95);
xlabel('FPSR (%)')
ylabel('GO categories reported as significant (%)')
legend(legendEntries,'Location','northwest')
title(whatSpecies)
colors = BF_getcmap('browngreen',3,1);
bh(1).FaceColor = colors{1};
bh(2).FaceColor = colors{3};
ax = gca;
ax.XLim = [0,30];
f.Position = [1058         972         460         219];

% bar(binMeans/10000,FPSR_repOrNot,'stacked')


%-------------------------------------------------------------------------------
