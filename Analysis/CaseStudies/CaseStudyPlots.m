% First filter on results we're interested in plotting:
justCortex = true;
showZero = false;
sigThresh = 0.05;

%-------------------------------------------------------------------------------
% Load and process data:
%-------------------------------------------------------------------------------
fprintf(1,'Loading case study results to file\n');
load('AllCaseStudies.mat','resultsTables')
[resultsTables,numSig] = ProcessCaseStudy(resultsTables,justCortex,sigThresh,false);
theAnalyses = fieldnames(resultsTables);
numAnalyses = length(theAnalyses);
fprintf(1,'%u enrichment analyses (across nulls) to analyze\n',numAnalyses);

%===============================================================================
% Show the number of significant hits as a bar chart
%===============================================================================
if ~showZero
    myFilter = sum(numSig,2)>0;
    numSigShow = numSig(myFilter,:);
    theAnalysesShow = theAnalyses(myFilter);
end
f = figure('color','w');
ax = gca();
[~,theSort] = sort(numSigShow(:,1),'descend');
bar(numSigShow(theSort,:))
ax.XTick = 1:size(numSigShow,1);
ax.XTickLabel = theAnalysesShow(theSort);
ax.XTickLabelRotation = 30;
ax.TickLabelInterpreter = 'None';
legend('random-Gene','SBP-rand','SBP-spatial')
title('Significant enrichment across phenotypes')
xlabel('Enrichment analysis')
ylabel('Number of significant hits (FDR_N < 0.05)')

%-------------------------------------------------------------------------------
% Zoom in on an interesting analysis?
%-------------------------------------------------------------------------------
% Show bar chart of some of the top hits: (bar for each null for some of the top
% significant categories)


%-------------------------------------------------------------------------------
% What about getting some picture of commonly significant categories?
%-------------------------------------------------------------------------------
% Retrieve the significant GO IDs from each analysis
nullsCorr = {'pValZCorrRandomGene','pValZCorrRandomMap','pValZCorrSpatialLag'};
numNulls = length(nullsCorr);

sigIDs = cell(numAnalyses,numNulls);
for i = 1:numAnalyses
    for k = 1:numNulls
        isSignificant = (resultsTables.(theAnalyses{i}).(nullsCorr{k}) < sigThresh);
        sigIDs{i,k} = resultsTables.(theAnalyses{i}).GOID(isSignificant);
    end
end

% Agglomerate:
allIDs_randomGene = vertcat(sigIDs{:,1});
allUniqueIDs = unique(allIDs_randomGene);
numOccurrence = arrayfun(@(x)sum(allIDs_randomGene==x),allUniqueIDs);
[~,ix] = max(numOccurrence);
theIDToPlot = allUniqueIDs(ix);

% Bar chart of one of these top categories:
f = figure('color','w');
