% EnrichmentCompareProcessingMethods()
% Idea is to compare different preprocessing steps on the enrichment of a given
% edge property
whatEdgeProperty = 'communicability';

%-------------------------------------------------------------------------------
% Fixed parameters:
energyOrDensity = 'energy'; % what gene expression data to use
pValOrStat = 'stat'; % 'pval','stat'
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
numIterationsErmineJ = 10000; % number of iterations for GSR in ermineJ

%-------------------------------------------------------------------------------
% Set up a structure array containing all of the different processing options:
absTypes = {false}; % false -> coexpression contribution increases with the statistic
corrTypes = {'Spearman'}; % {'Spearman','Pearson'};
normalizationGeneTypes = {'none','robustSigmoid'};
normalizationRegionTypes = {'none','zscore'};
correctDistanceTypes = {false};
pThresholds = [0.05];

cntr = 0;
for i = 1:length(absTypes)
    absType = absTypes{i};
    for j = 1:length(corrTypes)
        corrType = corrTypes{j};
        for k = 1:length(normalizationGeneTypes)
            normalizationGeneType = normalizationGeneTypes{k};
            for l = 1:length(normalizationRegionTypes)
                normalizationRegionType = normalizationRegionTypes{l};
                for m = 1:length(correctDistanceTypes)
                    correctDistance = correctDistanceTypes{m};
                    for n = 1:length(pThresholds)
                        pThreshold = pThresholds(n);
                        %----------------------------------------------------------
                        cntr = cntr + 1;
                        processingSteps(cntr).abs = absType;
                        processingSteps(cntr).corrType = corrType;
                        processingSteps(cntr).normalizationGene = normalizationGeneType;
                        processingSteps(cntr).normalizationRegion = normalizationRegionType;
                        processingSteps(cntr).correctDistance = correctDistance;
                        processingSteps(cntr).pThreshold = pThreshold;
                    end
                end
            end
        end
    end
end
numProcessingTypes = length(processingSteps);
theSwitches = fields(processingSteps(1));
numSwitches = length(theSwitches);
fprintf(1,'Comparing %u different processing parameters\n',numProcessingTypes);

%-------------------------------------------------------------------------------
% Get edge data:
%-------------------------------------------------------------------------------
C = load('Mouse_Connectivity_Data.mat'); % C stands for connectome data
f = figure('color','w');
for p = 1:length(pThresholds)
    A_bin = GiveMeAdj(C,'binary','ipsi',0,pThreshold);
    switch whatEdgeProperty
    case 'communicability'
        edgeData = communicability(A_bin);
        edgeData(~A_bin) = 0; % only put on real edges
    case 'betweenness'
        edgeData = edge_betweenness_bin(A_bin);
    case 'distance'
        edgeData = C.Dist_Matrix{1,1}/1000; % ipsilateral distances in the right hemisphere
    end
    subplot(2,length(pThresholds),2*(p-1)+1);
    histogram(edgeData(edgeData~=0));
    xlabel(whatEdgeProperty)
    d = C.Dist_Matrix{1,1}/1000; % ipsilateral distances in the right hemisphere
    isUpper = triu(true(size(d)),1);
    subplot(2,length(pThresholds),2*p);
    BF_PlotQuantiles(d(isUpper),edgeData(isUpper),15,false,false);
    % plot(d(isUpper),log10(edgeData(isUpper)),'.k');
    xlabel('d');
    ylabel(whatEdgeProperty)
end

%-------------------------------------------------------------------------------
% Get scores for genes:
%-------------------------------------------------------------------------------
enrichmentTables = cell(numProcessingTypes,1);
timer = tic;
for i = 1:numProcessingTypes
    fprintf(1,'%u/%u: norm-gene-%s, norm-reg-%s, corr-%s, dcorr-%u, pThresh-%.2f\n',i,numProcessingTypes,...
                                        processingSteps(i).normalizationGene,...
                                        processingSteps(i).normalizationRegion,...
                                        processingSteps(i).corrType,...
                                        processingSteps(i).correctDistance,...
                                        processingSteps(i).pThreshold);
    % Load in our gene data:
    [geneData,geneInfo,structInfo] = LoadMeG({processingSteps(i).normalizationGene,...
                        processingSteps(i).normalizationRegion},energyOrDensity);
    % Compute gene scores:
    % (sometimes entrez IDs change -- e.g., when matching to distance results)
    [gScore,geneEntrezIDs] = GiveMeGCC(edgeData,geneData,geneInfo.entrez_id,processingSteps(i).corrType,...
                    processingSteps(i).correctDistance,processingSteps(i).abs,thresholdGoodGene,pValOrStat);

    % Do enrichment:
    fileNameWrite = writeErmineJFile('tmp',gScore,geneEntrezIDs,whatEdgeProperty);
    ermineJResults = RunErmineJ(fileNameWrite,numIterationsErmineJ);
    % <<<Keep ermineJResults for p-values under 0.1>>>
    enrichmentTables{i} = ermineJResults(ermineJResults.pVal_corr < 0.1,:);
    % Give user feedback
    fprintf(1,'\n\n----%u/%u (%s remaining)\n\n',i,numProcessingTypes,...
                            BF_thetime((numProcessingTypes-i)*(toc(timer)/i)));
end


%-------------------------------------------------------------------------------
% Ok, and now we want to summarize the results
%-------------------------------------------------------------------------------

% 1. What are the different GO IDs implicated
GOIDs = cellfun(@(x)x.GOID,enrichmentTables,'UniformOutput',0);
allGOIDs = unique(vertcat(GOIDs{:}));
numGOIDs = length(allGOIDs);
% map to names
allGONames = cell(numGOIDs,1);
for i = 1:numGOIDs
    isHere = cellfun(@(x)ismember(allGOIDs{i},x),GOIDs);
    isHere = find(isHere,1,'first');
    thisRow = strcmp(enrichmentTables{isHere}.GOID,allGOIDs{i});
    allGONames{i} = enrichmentTables{isHere}.GOName{thisRow};
end
allGOLabels = arrayfun(@(x)sprintf('%s (%s)',allGONames{x},allGOIDs{x}),...
                            1:numGOIDs,'UniformOutput',false);

% 2. Prepare output
summaryTable = ones(numGOIDs,numProcessingTypes)*NaN;
for i = 1:numProcessingTypes
    if isempty(enrichmentTables{i}), continue; end
    [~,ia,ib] = intersect(allGOIDs,enrichmentTables{i}.GOID);
    summaryTable(ia,i) = enrichmentTables{i}.pVal(ib);
end
% order GO categories by relevance:
summaryTableSat = summaryTable;
summaryTableSat(isnan(summaryTable)) = 0.1;
propSig = mean(summaryTableSat,2);
[~,ix_GO] = sort(propSig,'descend');
propSig = mean(summaryTableSat,1);
[~,ix_pMeth] = sort(propSig,'descend');

labelTable = zeros(numSwitches,numProcessingTypes);
groupNamesAll = cell(numSwitches,1);
gidsAll = cell(numSwitches,1);
for i = 1:numSwitches
    if ischar(processingSteps(1).(theSwitches{i}))
        [gidsAll{i},groupNamesAll{i}] = findgroups({processingSteps.(theSwitches{i})});
    else
        [gidsAll{i},groupNamesAll{i}] = findgroups([processingSteps.(theSwitches{i})]);
    end
    labelTable(i,:) = gidsAll{i};
end

%-------------------------------------------------------------------------------
% Plot a table summarizing enrichment results
%-------------------------------------------------------------------------------

f = figure('color','w'); hold on;
maxSummaryTable = max(summaryTable(:));
subplot(5,1,1); ax1 = gca;
imagesc((BF_NormalizeMatrix(labelTable(:,ix_pMeth)','maxmin')'*0.999+1.001)*maxSummaryTable)
ax1.YLim = [0.5,0.5+numSwitches];
ax1.YTick = 1:numSwitches;
ax1.YTickLabel = theSwitches;
caxis([0,maxSummaryTable*2])
title(whatEdgeProperty)
subplot(5,1,2:5); hold on; ax2 = gca;
BF_imagesc(summaryTable(ix_GO,ix_pMeth));
ax2.XLim = [0.5,numProcessingTypes+0.5];
ax2.YLim = [0.5,numGOIDs + 0.5];
ax2.YTick = 1:numGOIDs;
ax2.YTickLabel = allGOLabels(ix_GO);
% plot([0.5,numProcessingTypes+0.5],(numGOIDs+0.5)*ones(2,1),'k')
ax2.XTick = 1:numProcessingTypes;
ax2.XTickLabel = [];
theMap = [BF_getcmap('blues',9,0);BF_getcmap('spectral',11,0)];
colormap(theMap(1:end-2,:))
caxis([0,maxSummaryTable*2])
% Label group names
subplot(5,1,1)
for i = 1:numSwitches
    for j = 1:numProcessingTypes
        if iscell(groupNamesAll{i})
            text(j,i,groupNamesAll{i}{gidsAll{i}(ix_pMeth(j))},'color','w')
        else
            text(j,i,num2str(groupNamesAll{i}(gidsAll{i}(ix_pMeth(j)))),'color','w')
        end
    end
    % numGroups = length(groupNamesAll{i});
    % for j = 1:numGroups
    %     if iscell(groupNamesAll{i})
    %         text(find(gidsAll{i}(ix_pMeth)==j,1,'first'),numGOIDs+i,groupNamesAll{i}{j},'color','w')
    %     else
    %         text(find(gidsAll{i}(ix_pMeth)==j,1,'first'),numGOIDs+i,num2str(groupNamesAll{i}(j)),'color','w')
    %     end
    % end
end
linkaxes([ax1,ax2],'x')
ax1.XLim = [0.5,numProcessingTypes+0.5];
ax2.Position = [0.4141    0.1100    0.5686    0.6423];
ax1.Position = [0.4141    0.7526    0.5686    0.1724];
