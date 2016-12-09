% Given a set of processing parameters, compares enrichment across
% different edge-based measures
% We'll use EnrichmentCompareProcessingMethods as a template
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Fixed processing parameters:
% Connectome processing
connectomeSource = 'Oh'; % 'Oh-cortex'
pThreshold = 0.05;
whatHemispheres = 'right';
justCortex = false;
% Gene processing
energyOrDensity = 'energy'; % what gene expression data to use
normalizationGene = 'none';
normalizationRegion = 'none'; %,'zscore'; % {'none','zscore'}
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
[A_bin,regionStruct] = GiveMeAdj(connectomeSource,pThreshold,true,whatHemispheres,justCortex);
A_wei = GiveMeAdj(connectomeSource,pThreshold,false,whatHemispheres,justCortex);
edgeMeasures = GiveMeAllEdgeMeasures(A_bin,A_wei);
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
    % <<<Keep ermineJResults for p-values under 0.1>>>
    enrichmentTables{i} = ermineJResults(ermineJResults.pVal_corr < enrichmentSigThresh,:);
    % Give user feedback
    fprintf(1,'\n\n----%u/%u (%s remaining)\n\n',i,numEdgeMeasures,...
                            BF_thetime((numEdgeMeasures-i)*(toc(timer)/i)));
end

%-------------------------------------------------------------------------------
% Prepare a table of the results
%-------------------------------------------------------------------------------
summaryTable = nan(numGOIDs,numProcessingTypes);
for i = 1:numProcessingTypes
    if isempty(enrichmentTables{i}), continue; end
    [~,ia,ib] = intersect(allGOIDs,enrichmentTables{i}.GOID);
    summaryTable(ia,i) = enrichmentTables{i}.pVal(ib);
end
% Order GO categories by relevance:
summaryTableSat = summaryTable;
summaryTableSat(isnan(summaryTable)) = 0.1;
propSig = mean(summaryTableSat,2);
[~,ix_GO] = sort(propSig,'descend');
propSig = mean(summaryTableSat,1);
[~,ix_pMeth] = sort(propSig,'descend');
%-------------------------------------------------------------------------------
% Plot a table summarizing enrichment results
%-------------------------------------------------------------------------------

f = figure('color','w'); hold on;
title(sprintf('(p_{FDR} < %.2f)',enrichmentSigThresh))
BF_imagesc(summaryTable(ix_GO,ix_pMeth));
ax2.XLim = [0.5,numProcessingTypes+0.5];
ax2.YLim = [0.5,numGOIDs + 0.5];
ax2.YTick = 1:numGOIDs;
ax2.YTickLabel = allGOLabels(ix_GO);
% plot([0.5,numProcessingTypes+0.5],(numGOIDs+0.5)*ones(2,1),'k')
ax2.XTick = 1:numProcessingTypes;
ax2.XTickLabel = [];
theMap = BF_getcmap('blues',9,0);
colormap(theMap(1:end-2,:))

linkaxes([ax1,ax2],'x')
ax1.XLim = [0.5,numProcessingTypes+0.5];
ax2.Position = [0.4141    0.1100    0.5686    0.6423];
ax1.Position = [0.4141    0.7526    0.5686    0.1724];
