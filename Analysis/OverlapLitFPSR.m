

% Load literature annotations:
load('LiteratureEnrichmentLoaded.mat','resultsTables','mouseOrHuman');

% Thresholds so that (i) only annotations with below-threshold significance are included, and
% (ii) only studies containing such below-threshold categories are filtered
% -> resultsTablesTh
pValCorrThreshold = 0.05;
resultsTablesTh = resultsTables;
resultsTablesTh = structfun(@(x)x(x.pValCorr<=pValCorrThreshold,:),resultsTablesTh,'UniformOutput',false);
emptyIndex = find(structfun(@(x)(height(x)==0),resultsTablesTh));
numEmpty = length(emptyIndex);
allFields = fieldnames(resultsTablesTh);
emptyFields = allFields(emptyIndex);
for i = 1:numEmpty
    resultsTablesTh = rmfield(resultsTablesTh,emptyFields{i});
end

% Convert to a table of GO categories, annotated by a list of studies
allFlaggedGOIDs = structfun(@(x)x.GOID,resultsTablesTh,'UniformOutput',false);
allFlaggedGOIDs = struct2cell(allFlaggedGOIDs);
allFlaggedGOIDs = unique(vertcat(allFlaggedGOIDs{:}));

% Now form a table:
% LiteratureGOTable = table();
numGOIDs = length(allFlaggedGOIDs);
studyList = cell(numGOIDs,1);
allFields = fieldnames(resultsTablesTh);
for i = 1:numGOIDs
    isListed = structfun(@(x)ismember(allFlaggedGOIDs(i),x.GOID),resultsTablesTh);
    studyList{i} = allFields(isListed);
end
GOID = allFlaggedGOIDs;
LiteratureGOTable = table(GOID,studyList);
numStudies = cellfun(@length,studyList);
[~,ix] = sort(numStudies,'descend');
LiteratureGOTable = LiteratureGOTable(ix,:);

%-------------------------------------------------------------------------------
% Now we'll get the mouse FPSR data:
whatSpecies = 'mouse';
numNullSamples_surrogate = 10000;
GOTable_FPSR = SurrogateEnrichmentProcess(whatSpecies,numNullSamples_surrogate,'randomUniform','');

% If at least X studies have quoted it at pCorr < 0.05, then we label it as a "reported category":
GOID_reported_1 = GOID(numStudies==1);
GOID_reported_2 = GOID(numStudies==2);
GOID_reported_3 = GOID(numStudies>2);
GOID_reported2More = GOID(numStudies>=2);

% Split into reported and unreported categories:
GOTable_notReported = GOTable_FPSR(~ismember(GOTable_FPSR.GOID,GOID),:);
GOTable_reported1 = GOTable_FPSR(ismember(GOTable_FPSR.GOID,GOID_reported_1),:);
GOTable_reported2 = GOTable_FPSR(ismember(GOTable_FPSR.GOID,GOID_reported_2),:);
GOTable_reported3 = GOTable_FPSR(ismember(GOTable_FPSR.GOID,GOID_reported_3),:);

GOTable_reported2More = GOTable_FPSR(ismember(GOTable_FPSR.GOID,GOID_reported2More),:);
GOTable_notreported2More = GOTable_FPSR(~ismember(GOTable_FPSR.GOID,GOID_reported2More),:);

% Distributions:
FPSR_reported = GOTable_reported2More.sumUnderSig;
FPSR_notReported = GOTable_notreported2More.sumUnderSig;
BF_JitteredParallelScatter({FPSR_reported,FPSR_notReported},true,true,true);

% Propotion by bins:
numBins = 10;
binEdges = linspace(0,max(GOTable_FPSR.sumUnderSig)+eps,numBins+1);
FPSR_repOrNot = zeros(numBins,4);
isInBin = @(x,Gtable) sum(Gtable.sumUnderSig>=binEdges(x) & Gtable.sumUnderSig < binEdges(x+1));
for i = 1:numBins
    FPSR_repOrNot(i,1) = isInBin(i,GOTable_notReported);
    FPSR_repOrNot(i,2) = isInBin(i,GOTable_reported1);
    FPSR_repOrNot(i,3) = isInBin(i,GOTable_reported2);
    FPSR_repOrNot(i,4) = isInBin(i,GOTable_reported3);
end
propReportedNorm = FPSR_repOrNot./sum(FPSR_repOrNot,2);
binMeans = mean([binEdges(1:end-1);binEdges(2:end)])
f = figure('color','w');
bar(100*binMeans/10000,100*propReportedNorm(:,3:end),'stacked')
xlabel('FPSR (%)')
ylabel('GO categories reported as significant (%)')
legend('2 studies','3+ studies')

% bar(binMeans/10000,FPSR_repOrNot,'stacked')


%-------------------------------------------------------------------------------

theStudyNames = fieldnames(resultsTables);
numEnrichmentStudies = length(theStudyNames);
for i = 1:numEnrichmentStudies
    % Look for the target in this study:
    if ismember(whatGOID,resultsTables.(theStudyNames{i}).GOID)
        fprintf(1,'Match for GOID %u in %s:     ',whatGOID,theStudyNames{i});
        here = (resultsTables.(theStudyNames{i}).GOID==whatGOID);
        fprintf(1,'pVal = %g\n',resultsTables.(theStudyNames{i}).pValCorr(here));
    end
end

end
