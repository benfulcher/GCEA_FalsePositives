% EnrichmentCompareProcessingMethods()
% Idea is to compare different preprocessing steps on the enrichment of a given
% edge property
whatEdgeProperty = 'betweenness';


%-------------------------------------------------------------------------------
% Fixed parameters:
energyOrDensity = 'energy'; % what gene expression data to use
pValOrStat = 'stat'; % 'pval','stat'
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
numIterationsErmineJ = 20000; % number of iterations for GSR in ermineJ

%-------------------------------------------------------------------------------
% Set up a structure array containing all of the different processing options:
absTypes = {true,false};
corrTypes = {'Pearson'}; % {'Spearman','Pearson'};
normalizationGeneTypes = {'none'}; % {'none','log10','robustSigmoid'};
normalizationRegionTypes = {'zscore'}; % {'none','zscore','robustSigmoid'};
correctDistanceTypes = {true,false};

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
                    %----------------------------------------------------------
                    cntr = cntr + 1;
                    processingSteps(cntr).abs = absType;
                    processingSteps(cntr).corrType = corrType;
                    processingSteps(cntr).normalizationGene = normalizationGeneType;
                    processingSteps(cntr).normalizationRegion = normalizationRegionType;
                    processingSteps(cntr).correctDistance = correctDistance;
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
pThreshold = 0.05;
A_bin = GiveMeAdj(C,'binary','ipsi',0,pThreshold);
edgeData = edge_betweenness_bin(A_bin);
f = figure('color','w');
histogram(edgeData(edgeData~=0));
xlabel('Edge betweenness (binary)')

%-------------------------------------------------------------------------------
% Get scores for genes:
%-------------------------------------------------------------------------------
enrichmentTables = cell(numProcessingTypes,1);
timer = tic;
for i = 1:numProcessingTypes
    % Load in our gene data:
    [GeneStruct,geneData] = LoadMeG(true,{processingSteps(i).normalizationGene,...
                            processingSteps(i).normalizationRegion},energyOrDensity);
    geneEntrezIDs = [GeneStruct.gene_entrez_id];
    % Compute gene scores:
    % (sometimes entrez IDs change -- e.g., when matching to distance results)
    [gScore,geneEntrezIDs] = GiveMeGCC(edgeData,geneData,geneEntrezIDs,processingSteps(i).corrType,...
                    processingSteps(i).correctDistance,thresholdGoodGene,pValOrStat);
    % Do enrichment:
    fileNameWrite = writeErmineJFile('tmp',gScore,geneEntrezIDs,'bin_betw');
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
f = figure('color','w'); hold on; ax = gca;
maxSummaryTable = max(summaryTable(:));
BF_imagesc([summaryTable;labelTable*maxSummaryTable+1e-6]);
ax.YLim = [0.5,numGOIDs + numSwitches+0.5];
ax.YTick = 1:numGOIDs + numSwitches;
ax.YTickLabel(1:numGOIDs) = allGOLabels;
ax.YTickLabel(numGOIDs+1:numGOIDs+numSwitches) = theSwitches;
plot([0.5,numProcessingTypes+0.5],(numGOIDs+0.5)*ones(2,1),'k')
ax.XTick = 1:numProcessingTypes;
colormap([BF_getcmap('blues',9,0);BF_getcmap('spectral',9,0)])
% Label group names
for i = 1:numSwitches
    numGroups = length(groupNamesAll{i});
    for j = 1:numGroups
        if iscell(groupNamesAll{i})
            text(find(gidsAll{i}==j,1,'first'),numGOIDs+i,groupNamesAll{i}{j},'color','w')
        else
            text(find(gidsAll{i}==j,1,'first'),numGOIDs+i,num2str(groupNamesAll{i}(j)),'color','w')
        end
    end
end
colorbar()
