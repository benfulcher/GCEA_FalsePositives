%-------------------------------------------------------------------------------
% Aim is to try to score gene categories across a range of null networks
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Fixed processing parameters:
%-------------------------------------------------------------------------------
% Connectome processing
connectomeSource = 'Oh';
pThreshold = 0.05;
whatHemispheres = 'right';
justCortex = false;
whatEdgeMeasure = 'ktotktot';
onlyOnEdges = true; % whether to put values only on existing edges
                    % (rather than all node pairs for some measures)
numNulls = 20;

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


%-------------------------------------------------------------------------------
% Compute randomized
fprintf(1,'Randomizing edges uniformly\n');
N = length(A_bin);
numLinks = sum(A_bin(:));

edgeMeasures = zeros(N,numNulls);
for i = 1:numNulls
    A_rand_vector = zeros(N*(N-1),1);
    rp = randperm(N*(N-1));
    A_rand_vector(rp(1:numLinks)) = 1;
    A_bin_i = squareform(A_rand_vector(1:end/2));
    A_bin_ii = squareform(A_rand_vector(end/2+1:end));
    upper = triu(true(size(A_bin)),+1);
    A_rand = zeros(size(A_bin));
    A_rand(upper) = A_bin_i(upper);
    lower = tril(true(size(A_bin)),-1);
    A_rand(lower) = A_bin_ii(lower);

    weightVector = A_wei(A_wei > 0);
    A_wei_rand = A_rand;
    A_wei_rand(A_wei_rand>0) = weightVector;

    ktot = sum(A_rand,1)' + sum(A_rand,2);
    product = ktot*ktot';
    product(A_rand == 0) = 0;
    edgeMeasures(:,i) = product;
end


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
% Don't regress distance:
distanceRegressor = [];

%-------------------------------------------------------------------------------
% Get scores for genes:
%-------------------------------------------------------------------------------
enrichmentTables = cell(numEdgeMeasures,1);
enrichmentSigThresh = 0.05;
fprintf(1,'---Interested in Biological Processes with FDR p < %g\n',enrichmentSigThresh);
timer = tic;
for i = 1:numNulls
    fprintf(1,'%u/%u: %s\n\n',i,numNulls,edgeMeasureNames{i});

    % Compute gene scores:
    % (sometimes entrez IDs change -- e.g., when matching to distance results)
    [gScore,geneEntrezIDs] = GiveMeGCC(edgeMeasures(:,i),geneData,...
                                geneInfo.entrez_id,corrType,distanceRegressor,absType,...
                                thresholdGoodGene,pValOrStat);

    % Do enrichment:
    

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
if numGOIDs > 0
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
end

%-------------------------------------------------------------------------------
% Look into the correlations within each GO group
