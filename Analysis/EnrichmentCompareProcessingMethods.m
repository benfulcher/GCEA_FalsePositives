% EnrichmentCompareProcessingMethods()
% Idea is to compare different preprocessing steps on the enrichment of a given
% edge property
whatEdgeProperty = 'bin-communicability';

%-------------------------------------------------------------------------------
% Fixed parameters:
energyOrDensity = 'energy'; % what gene expression data to use
pValOrStat = 'stat'; % 'pval','stat'
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
numIterationsErmineJ = 20000; % number of iterations for GSR in ermineJ
numQuantiles = 15; % quantiles with which to learn the distance relationship

%-------------------------------------------------------------------------------
% Set up a structure array containing all of the different processing options:
connectomeTypes = {'Oh-brain'}; % 'Oh-cortex'
absTypes = {'pos'}; % 'pos','neg','abs' -> e.g., pos -> coexpression contribution increases with the statistic
corrTypes = {'Spearman'}; % {'Spearman','Pearson'};
normalizationGeneTypes = {'none'};
normalizationRegionTypes = {'none','zscore'}; % {'none','zscore'}
correctDistanceTypes = {false,true};
pThresholds = [0.05];
processingSteps = struct();

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
                        for o = 1:length(connectomeTypes)
                            connectomeType = connectomeTypes{o};
                            %----------------------------------------------------------
                            cntr = cntr + 1;
                            processingSteps(cntr).absType = absType;
                            processingSteps(cntr).corrType = corrType;
                            processingSteps(cntr).normalizationGene = normalizationGeneType;
                            processingSteps(cntr).normalizationRegion = normalizationRegionType;
                            processingSteps(cntr).correctDistance = correctDistance;
                            processingSteps(cntr).pThreshold = pThreshold;
                            processingSteps(cntr).connectomeType = connectomeType;
                        end
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
% ~~~~~Plotting now -> DistanceDependenceEdge.m~~~~~~~~~~
C = load('Mouse_Connectivity_Data.mat','Dist_Matrix');
f = figure('color','w');
edgeData = cell(length(pThresholds),2);
for p = 1:length(pThresholds)
    switch connectomeTypes{1}
    case 'Oh-brain'
        A_bin = GiveMeAdj('Oh',pThresholds(p),true,'right',false);
        A_wei = GiveMeAdj('Oh',pThresholds(p),false,'right',false);
    case 'Oh-cortex'
        A_bin = GiveMeAdj('Oh',pThresholds(p),true,'right',true);
        A_wei = GiveMeAdj('Oh',pThresholds(p),false,'right',true);
    otherwise
        error('Unknown connectome: %s',connectomeTypes);
    end

    % Compute the edge data:
    onlyOnEdges = true;
    if strcmp(whatEdgeProperty,'distance')
        edgeData{p,1} = C.Dist_Matrix{1,1}/1000; % ipsilateral distances in the right hemisphere
    else
        edgeData{p,1} = GiveMeEdgeMeasure(whatEdgeProperty,A_bin,A_wei,onlyOnEdges);
    end

    %---------------------------------------------------------------------------
    subplot(2,length(pThresholds),2*(p-1)+1);
    histogram(edgeData{p,1}(edgeData{p,1}~=0));
    xlabel(whatEdgeProperty)
    d = C.Dist_Matrix{1,1}/1000; % ipsilateral distances in the right hemisphere
    title(pThresholds(p))
    subplot(2,length(pThresholds),2*p);
    connValues = edgeData{p,1} > 0;
    edgeDataCorrected = BF_PlotQuantiles(d(connValues),edgeData{p,1}(connValues),numQuantiles,false,false);
    title(sprintf('%s on %u edges',whatEdgeProperty,sum(connValues(:))))
    xlabel('d');
    ylabel(whatEdgeProperty)
    edgeData{p,2} = zeros(size(edgeData{p,1}));
    edgeData{p,2}(connValues) = edgeDataCorrected;
end
drawnow

%-------------------------------------------------------------------------------
% Get scores for genes:
%-------------------------------------------------------------------------------
enrichmentTables = cell(numProcessingTypes,1);
enrichmentSigThresh = 0.05;
fprintf(1,'---Interested in Biological Processes with FDR p < %g\n',enrichmentSigThresh);
timer = tic;
for i = 1:numProcessingTypes
    fprintf(1,'%u/%u: %s: norm-gene-%s, norm-reg-%s, corr-%s, dcorr-%u, pThresh-%.2f, abs-%s\n',...
                                        i,numProcessingTypes,...
                                        processingSteps(i).connectomeType,...
                                        processingSteps(i).normalizationGene,...
                                        processingSteps(i).normalizationRegion,...
                                        processingSteps(i).corrType,...
                                        processingSteps(i).correctDistance,...
                                        processingSteps(i).pThreshold,...
                                        processingSteps(i).absType);

    % Compute edge-level statistics:
    % [a quantile-based distance correction could appear here [false]; but now
    % we're doing a partial correlation-based correction]
    [edgeData,regionStruct] = GiveMeEdgeStat(processingSteps(i).connectomeType,processingSteps(i).pThreshold,...
                            whatEdgeProperty,false,numQuantiles);

    % Load in our gene data, properly processed:
    [geneData,geneInfo,structInfo] = LoadMeG({processingSteps(i).normalizationGene,...
                        processingSteps(i).normalizationRegion},energyOrDensity);

    % Check gene data matches connectome data
    if ~all([regionStruct.id]'==structInfo.id)
        % Take subset
        [~,ia,ib] = intersect([regionStruct.id]',structInfo.id,'stable');
        geneData = geneData(ib,:);
        structInfo = structInfo(ib,:);
        fprintf(1,'Gene data matched to subset of %u Allen regions\n',length(ib));
    end

    % Correct distance?
    if processingSteps(i).correctDistance
        C = load('Mouse_Connectivity_Data.mat','Dist_Matrix');
        distanceRegressor = C.Dist_Matrix{1,1};
        fprintf(1,'Regressing ipsilateral distances\n');
    else
        distanceRegressor = []; % just compute normal correlations
    end

    % Compute gene scores:
    % (sometimes entrez IDs change -- e.g., when matching to distance results)
    [gScore,geneEntrezIDs] = GiveMeGCC(edgeData,geneData,geneInfo.entrez_id,...
                                    processingSteps(i).corrType,distanceRegressor,...
                                    processingSteps(i).absType,thresholdGoodGene,pValOrStat);

    % Do enrichment:
    fileNameWrite = writeErmineJFile('tmp',gScore,geneEntrezIDs,whatEdgeProperty);
    ermineJResults = RunErmineJ(fileNameWrite,numIterationsErmineJ);
    % <<<Keep ermineJResults for p-values under 0.1>>>
    enrichmentTables{i} = ermineJResults(ermineJResults.pVal_corr < enrichmentSigThresh,:);
    % Give user feedback
    fprintf(1,'\n\n----%u/%u (%s remaining)\n\n',i,numProcessingTypes,...
                            BF_thetime((numProcessingTypes-i)*(toc(timer)/i)));
end


%-------------------------------------------------------------------------------
% Ok, and now we want to summarize the results
%-------------------------------------------------------------------------------
doReorder = true;
[summaryTable,allGOLabels,allGONames,allGOIDs,ix_pMeth] = PrepareSummaryTable(enrichmentTables,doReorder);


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
title(sprintf('%s (p_{FDR} < %.2f)',whatEdgeProperty,enrichmentSigThresh))
subplot(5,1,2:5); hold on; ax2 = gca;
BF_imagesc(summaryTable(:,ix_pMeth));
ax2.XLim = [0.5,numProcessingTypes+0.5];
ax2.YLim = [0.5,numGOIDs + 0.5];
ax2.YTick = 1:numGOIDs;
ax2.YTickLabel = allGOLabels(:);
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
