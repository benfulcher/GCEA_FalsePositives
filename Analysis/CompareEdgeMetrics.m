% Given a set of processing parameters, compares enrichment across
% different edge-based measures
% We'll use EnrichmentCompareProcessingMethods as a template
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Fixed processing parameters:
%-------------------------------------------------------------------------------
% Connectome processing
connectomeSource = 'Oh'; % 'Oh-cortex'
pThreshold = 0.05;
whatHemispheres = 'right';
justCortex = false;
onlyOnEdges = true; % whether to put values only on existing edges
                    % (rather than all node pairs for some measures)
randomizeEdges = true; % whether to randomize edges [doesn't work yet]

% Gene processing
energyOrDensity = 'energy'; % what gene expression data to use
normalizationGene = 'none'; % 'none', 'mixedSigmoid'
normalizationRegion = 'none'; % 'none', 'zscore'

% Correlations:
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
corrType = 'Spearman'; % {'Spearman','Pearson'};
pValOrStat = 'stat'; % 'pval','stat'
absType = 'pos'; % 'pos','neg','abs' -> e.g., pos -> coexpression contribution increases with the statistic
correctDistance = false; % false,true;

% Enrichment settings:
numIterationsErmineJ = 20000; % number of iterations for GSR in ermineJ

%-------------------------------------------------------------------------------

% Define a set of edge measures to compare:
[A_bin,regionStruct,adjPVals] = GiveMeAdj(connectomeSource,pThreshold,true,whatHemispheres,justCortex);
A_wei = GiveMeAdj(connectomeSource,pThreshold,false,whatHemispheres,justCortex);
edgeMeasures = GiveMeAllEdgeMeasures(A_bin,A_wei,onlyOnEdges,adjPVals);
edgeMeasureNames = fieldnames(edgeMeasures);
numEdgeMeasures = length(edgeMeasureNames);
fprintf(1,'Comparing %u edge measures\n',numEdgeMeasures);

%-------------------------------------------------------------------------------
% Get, and match gene data
[geneData,geneInfo,structInfo] = LoadMeG({normalizationGene,normalizationRegion},energyOrDensity);
% Check gene data matches connectome data
if ~all([regionStruct.id]'==structInfo.id)
    % Take subset
    [~,ia,ib] = intersect([regionStruct.id]',structInfo.id,'stable');
    geneData = geneData(ib,:);
    structInfo = structInfo(ib,:);
    fprintf(1,'Gene data matched to subset of %u Allen regions\n',length(ib));
end

%-------------------------------------------------------------------------------
% Get distance to use for regressor?:
if correctDistance
    C = load('Mouse_Connectivity_Data.mat','Dist_Matrix');
    distanceRegressor = C.Dist_Matrix{1,1};
    fprintf(1,'Regressing ipsilateral distances\n');
else
    distanceRegressor = []; % just compute normal correlations
end

%-------------------------------------------------------------------------------
% Randomize to isolate the driving effect:
if randomizeEdges
    % fprintf(1,'RANDOMIZING THE ASSIGNMENT OF GENE DATA TO REGIONS!!\n');
    % ix = randperm(length(A_bin));
    % geneData = geneData(ix,:);
    fprintf(1,'Assinging gene information randomly!\n');
    ix = randperm(size(geneData,2));
    geneInfo = geneInfo(ix,:);
    % A_bin = A_bin(ix,ix);
    % A_wei = A_wei(ix,ix);
    % This randomizes the topology under connectivity constraints
    % numIter = 1000;
    % A_bin = randomizeTopology(A_bin,numIter);
    % A_wei = randomizeTopology(A_wei,numIter);
    % fprintf(1,' Done!!\n');
end

%-------------------------------------------------------------------------------
% Get scores for genes:
%-------------------------------------------------------------------------------
enrichmentTables = cell(numEdgeMeasures,1);
enrichmentSigThresh = 0.05;
fprintf(1,'---Interested in Biological Processes with FDR p < %g\n',enrichmentSigThresh);
timer = tic;
for i = 1:numEdgeMeasures
    fprintf(1,'%u/%u: %s\n\n',i,numEdgeMeasures,edgeMeasureNames{i});

    % Compute gene scores:
    % (sometimes entrez IDs change -- e.g., when matching to distance results)
    [gScore,geneEntrezIDs] = GiveMeGCC(edgeMeasures.(edgeMeasureNames{i}),geneData,...
                                geneInfo.entrez_id,corrType,distanceRegressor,absType,...
                                thresholdGoodGene,pValOrStat);

    % Do enrichment:
    fileNameWrite = writeErmineJFile('tmp',gScore,geneEntrezIDs,edgeMeasureNames{i});
    ermineJResults = RunErmineJ(fileNameWrite,numIterationsErmineJ);
    % <<<Keep ermineJResults for p-values under (e.g., 0.1)>>>
    enrichmentTables{i} = ermineJResults(ermineJResults.pVal_corr < enrichmentSigThresh,:);
    % Give user feedback
    fprintf(1,'\n\n----%u/%u (%s remaining)\n\n',i,numEdgeMeasures,...
                            BF_thetime((numEdgeMeasures-i)*(toc(timer)/i)));
end

%-------------------------------------------------------------------------------
% Save
%-------------------------------------------------------------------------------
labelName = sprintf('%s-%s-%s-%s-%s-dist%u.mat',connectomeSource,...
                    normalizationGene,normalizationRegion,corrType,absType,...
                    correctDistance);
matName = [labelName,'.mat'];
fprintf(1,'Saving all outputs to %s...',matName);
save(fullfile(pwd,'DataOutputs',matName));
fprintf(1,' Saved.\n');

%-------------------------------------------------------------------------------
% Prepare a table of the results
%-------------------------------------------------------------------------------
doReorder = true;
[summaryTable,allGOLabels,allGONames,allGOIDs,ix_edgeM] = PrepareSummaryTable(enrichmentTables,doReorder);
numGOIDs = length(allGOIDs);
%-------------------------------------------------------------------------------
% Plot a table summarizing enrichment results
%-------------------------------------------------------------------------------
f = figure('color','w'); hold on; ax = gca;
title(sprintf('%s (p_{FDR} < %.2f)',labelName,enrichmentSigThresh))
BF_imagesc(summaryTable);
ax.XLim = [0.5,numEdgeMeasures+0.5];
ax.YLim = [0.5,numGOIDs + 0.5];
ax.YTick = 1:numGOIDs;
ax.YTickLabel = allGOLabels;
ax.XTickLabelRotation = 90;
ax.TickLabelInterpreter = 'none';
% plot([0.5,numEdgeMeasures+0.5],(numGOIDs+0.5)*ones(2,1),'k')
ax.XTick = 1:numEdgeMeasures;
ax.XTickLabel = edgeMeasureNames(ix_edgeM);
colormap(BF_getcmap('blues',9,0));
