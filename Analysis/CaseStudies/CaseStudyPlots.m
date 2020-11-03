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

nullsRaw = {'pValZRandomGene','pValZRandomMap','pValZSpatialLag'};
nullsCorr = {'pValZCorrRandomGene','pValZCorrRandomMap','pValZCorrSpatialLag'};
numNulls = length(nullsCorr);

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
f.Position(3:4) = [623,376];
fileName = 'sigHitsPhenotypes.svg';
saveas(f,fileName,'svg')
fprintf(1,'Saved to %s\n',fileName);

%-------------------------------------------------------------------------------
% Zoom in on an interesting analysis?
%-------------------------------------------------------------------------------
% Show bar chart of some of the top hits: (bar for each null for some of the top
% significant categories)
% theAnalysisToPlot = 'mouseCortexdegree';
theAnalysisToPlot = 'mouseCortex_oligodendrocytes';

doCorrected = false;
theResultsTable = resultsTables.(theAnalysisToPlot);
sigCategoryIndices = find(theResultsTable.pValZCorrRandomGene < sigThresh);
numCategories = length(sigCategoryIndices);
theTable = nan(numCategories,numNulls);
theCategoryNames = cell(numCategories,1);
for a = 1:numCategories
    for n = 1:numNulls
        if doCorrected
            theTable(a,n) = theResultsTable.(nullsCorr{n})(sigCategoryIndices(a));
        else
            theTable(a,n) = theResultsTable.(nullsRaw{n})(sigCategoryIndices(a));
        end
        if n==1
            theCategoryNames{a} = theResultsTable.GOName{sigCategoryIndices(a)};
        end
    end
end

%-------------------------------------------------------------------------------
f = figure('color','w');
ax = gca();
hold('on')
bar(-log10(theTable),'Horizontal','on')
ax.YTick = 1:numCategories;
ax.YTickLabel = theCategoryNames;
ax.TickLabelInterpreter = 'None';
if doCorrected
    plot(ones(2,1)*-log10(sigThresh),ax.YLim,'--r');
else
    theApproxThresh = max(theResultsTable.(nullsRaw{1})(sigCategoryIndices));
    plot(ones(2,1)*-log10(theApproxThresh),ax.YLim,'--r');
end
xlabel('-log10(p)')
ylabel('Analysis')
title(theAnalysisToPlot,'interpreter','none')
f.Position = [1000,1070,663,268];

%-------------------------------------------------------------------------------
% What about getting some picture of commonly significant categories?
%-------------------------------------------------------------------------------
% Retrieve the significant GO IDs from each analysis
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
theTable = nan(numAnalyses,numNulls);
for a = 1:numAnalyses
    for n = 1:numNulls
        theIndex = resultsTables.(theAnalyses{a}).GOID==theIDToPlot;
        theTable(a,n) = resultsTables.(theAnalyses{a}).(nullsCorr{n})(theIndex);
        if a==1 && n==1
            theCategoryName = resultsTables.(theAnalyses{a}).GOName{theIndex};
        end
    end
end
% Filter by significant
toPlot = min(theTable,[],2) < sigThresh;
theTable = theTable(toPlot,:);
theAnalysesFilter = theAnalyses(toPlot);

f = figure('color','w');
ax = gca();
hold('on')
bar(-log10(theTable),'Horizontal','on')
ax.YTick = 1:length(theAnalysesFilter);
ax.YTickLabel = theAnalysesFilter;
ax.TickLabelInterpreter = 'None';
plot(ones(2,1)*-log10(sigThresh),ax.YLim,'--r');
xlabel('-log10(p)')
ylabel('Analysis')
title(theCategoryName)
