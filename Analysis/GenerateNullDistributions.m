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
whatEdgeWeight = 'NCD'; % 'NCD' (normalized connection density), 'CS' (connection strength)
whatEdgeMeasure = 'wei_communicability'; % 'bin_edgeBet', 'ktot_ktot', 'wei_communicability'
onlyOnEdges = true; % whether to put values only on existing edges
                    % (rather than all node pairs for some measures)
useFakeConnectome = false; % e.g., use a connectome with 50% link density

% Randomization
randomizeHow = 'permutedGeneDep'; % 'topology', 'uniformTopology', 'permutedGeneDep', 'shuffleEdgeVals'
numNulls = 100;

% Gene processing
energyOrDensity = 'energy'; % what gene expression data to use
normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
normalizationRegion = 'none'; % 'none', 'zscore'
subsetOfGenes = []; % only look at the first X genes. Set to empty for all
                     % genes & save a .mat file of results

% GO settings
processFilter = 'biological_process';
sizeFilter = [5,200];

% Correlations:
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
corrType = 'Spearman'; % 'Spearman','Pearson'
pValOrStat = 'stat'; % 'pval','stat'
absType = 'pos'; % 'pos','neg','abs' -> e.g., pos -> coexpression contribution increases with the statistic
correctDistance = true; % false,true;


%-------------------------------------------------------------------------------
% Define a set of edge measures to compare:
if ~useFakeConnectome
    [A_bin,regionStruct,adjPVals] = GiveMeAdj(connectomeSource,pThreshold,true,...
                                    whatEdgeWeight,whatHemispheres,justCortex);
    A_wei = GiveMeAdj(connectomeSource,pThreshold,false,...
                                    whatEdgeWeight,whatHemispheres,justCortex);
else
    % A matrix with 50% connection probability:
    numEdgesUpper = (213*212)/2;
    edgesUpperEqual = [ones(numEdgesUpper/2,1);zeros(numEdgesUpper/2,1)];
    A_bin = zeros(213);
    A_bin(triu(true(size(A_bin)),+1)) = edgesUpperEqual;
    A_wei = A_bin;
end

%-------------------------------------------------------------------------------
% Compute randomized
switch randomizeHow
case 'uniformTopology'
    % Randomize the topology to have a uniform probability across all possible edges
    fprintf(1,'Randomizing edges uniformly\n');
    N = length(A_bin);
    numLinks = sum(A_bin(:));

    edgeMeasures = cell(numNulls+1,1);
    for i = 1:numNulls+1
        if i == 1
            A_rand = A_bin;
        else
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
        end

        weightVector = A_wei(A_wei > 0);
        A_wei_rand = A_rand;
        A_wei_rand(A_wei_rand > 0) = weightVector;

        % Compute the desired edge measure:
        edgeMeasures{i} = GiveMeEdgeMeasure(whatEdgeMeasure,A_rand,A_wei_rand,onlyOnEdges);
    end

case 'topology'
    numIter = 50; % (num randomizations per edge)
    fprintf(1,'Generating %u degree (+in-strength) topological nulls using %u randomizations per edge\n',...
                        numNulls,numIter);
    edgeMeasures = cell(numNulls+1,1);
    for i = 1:numNulls+1
        if i == 1
            A_rand = A_bin;
            A_wei_rand = A_wei;
        else
            A_rand = randmio_dir(A_bin,numIter);
            A_wei_rand = randmio_dir(A_wei,numIter);
        end
        % Compute the desired edge measure:
        edgeMeasures{i} = GiveMeEdgeMeasure(whatEdgeMeasure,A_rand,A_wei_rand,onlyOnEdges);
    end

case 'permutedGeneDep'
    % Permute gene profiles assigned to regions (later)
    % Edge measure stays constant using the correct topology:
    edgeMeasures = GiveMeEdgeMeasure(whatEdgeMeasure,A_bin,A_wei,onlyOnEdges);

case 'shuffleEdgeVals'
    % Shuffle values given to each edge
    % First compute the proper values:
    edgeMeasure0 = GiveMeEdgeMeasure(whatEdgeMeasure,A_bin,A_wei,onlyOnEdges);
    % Get the vector of values:
    edgeValsVector = edgeMeasure0(edgeMeasure0 > 0);
    % Now shuffle the values many times across the edges, preserving the binary topology
    edgeMeasures = cell(numNulls+1,1);
    for i = 1:numNulls+1
        if i==1
            edgeMeasures{i} = edgeMeasure0;
        else
            edgeMeasures{i} = zeros(size(A_bin));
            edgeMeasures{i}(A_bin > 0) = edgeValsVector(randperm(length(edgeValsVector)));
        end
    end
end

%-------------------------------------------------------------------------------
% Retrieve and match gene data
[geneData,geneInfo,structInfo] = LoadMeG({normalizationGene,normalizationRegion},energyOrDensity);
% SUBSET FIRST X GENES:
if ~isempty(subsetOfGenes)
    warning('Only looking at a random set of %u genes',subsetOfGenes);
    rp = randperm(size(geneData,2));
    rp = rp(1:subsetOfGenes);
    geneData = geneData(:,rp);
    geneInfo = geneInfo(rp,:);
end
% Check gene data matches connectome data
if ~all([regionStruct.id]'==structInfo.id)
    % Take subset
    [~,ia,ib] = intersect([regionStruct.id]',structInfo.id,'stable');
    geneData = geneData(ib,:);
    structInfo = structInfo(ib,:);
    fprintf(1,'Gene data matched to subset of %u Allen regions\n',length(ib));
end
[numRegions,numGenes] = size(geneData);

% Regress distance:
if correctDistance
    C = load('Mouse_Connectivity_Data.mat','Dist_Matrix');
    distanceRegressor = C.Dist_Matrix{1,1};
else
    distanceRegressor = [];
end

%-------------------------------------------------------------------------------
% Get GO data
% (include only annotations for genes with entrez IDs that are in our dataset)
[GOTable,geneEntrezAnnotations] = GetFilteredGOData(processFilter,sizeFilter,geneInfo.entrez_id);
sizeGOCategories = cellfun(@length,geneEntrezAnnotations);
numGOCategories = height(GOTable);

%-------------------------------------------------------------------------------
% Get scores for genes:
%-------------------------------------------------------------------------------
categoryScores = nan(numGOCategories,numNulls+1);
gScores = cell(numNulls+1,1);
entrezIDsKept = cell(numNulls+1,1);
enrichmentSigThresh = 0.05;
fprintf(1,'---Interested in Biological Processes with FDR p < %g\n',enrichmentSigThresh);
timer = tic;
for i = 1:numNulls+1
    fprintf(1,'%u/%u\n\n',i,numNulls+1);

    % Compute gene scores:
    % (sometimes entrez IDs change -- e.g., when matching to distance results)
    switch randomizeHow
    case {'topology','uniformTopology','shuffleEdgeVals'}
        % each null has a different set of edge measures:
        theEdgeData = edgeMeasures{i};
        theGeneData = geneData;
    case 'permutedGeneDep'
        % each null should use the same edge measure but different permutations of the gene data
        theEdgeData = edgeMeasures;
        if i==1
            % (don't permute the first one!)
            rp = 1:numRegions;
        else
            rp = randperm(numRegions);
        end
        theGeneData = geneData(rp,:);
    end
    % Compute the score:
    if strcmp(whatEdgeMeasure,'connected')
        % We want to score each gene by how much (relative to expectation of distance)
        % it's GCC scores differ
        [gScores{i},entrezIDsKept{i}] = ConnectedGCC(theEdgeData,theGeneData,...
                            geneInfo.entrez_id,corrType,distanceRegressor,absType,...
                            thresholdGoodGene,pValOrStat);
    else
        [gScores{i},entrezIDsKept{i}] = GiveMeGCC(theEdgeData,theGeneData,...
                            geneInfo.entrez_id,corrType,distanceRegressor,absType,...
                            thresholdGoodGene,pValOrStat);
    end

    % Record mean scores for each category:
    for j = 1:numGOCategories
        matchMe = ismember(entrezIDsKept{i},geneEntrezAnnotations{j});
        if sum(matchMe) <= 1
            continue
        end
        categoryScores(j,i) = nanmean(gScores{i}(matchMe));
    end

    % Give user feedback
    fprintf(1,'\n\n----%u/%u (%s remaining)\n\n',i,numNulls+1,...
                            BF_thetime((numNulls+1-i)*(toc(timer)/i)));
end

%-------------------------------------------------------------------------------
% Compute p-values
%-------------------------------------------------------------------------------
fprintf(1,'GO categories with correlations to %s than %s nulls\n',whatEdgeMeasure,randomizeHow);
whatTail = 'right';
[meanNull,stdNull,pValsPerm,pValsZ,pValsZ_corr] = EstimatePVals(categoryScores,numNulls,whatTail);

%-------------------------------------------------------------------------------
% Save to mat file:
%-------------------------------------------------------------------------------
if isempty(subsetOfGenes)
    fileName = sprintf('%s-%s-%s-G%s_R%s-%unulls.mat',whatEdgeMeasure,randomizeHow,...
                    processFilter,normalizationGene,normalizationRegion,numNulls);
    save(fullfile('DataOutputs',fileName));
    fprintf(1,'Saved %s\n',fileName);
end

%-------------------------------------------------------------------------------
% List categories with lowest p-values, or highest mean across nulls, etc.
%-------------------------------------------------------------------------------
ListCategories(geneInfo,GOTable,geneEntrezAnnotations,meanNull,pValsZ,pValsZ_corr);

%-------------------------------------------------------------------------------
% Check that the mean null score for each gene is zero
%-------------------------------------------------------------------------------
f = figure('color','w');
allKeptEntrez = unique(vertcat(entrezIDsKept{:}));
gScoresMat = nan(length(allKeptEntrez),numNulls);
for i = 1:numNulls
    [~,~,perm] = intersect(allKeptEntrez,entrezIDsKept{i+1},'stable');
    gScoresMat(:,i) = gScores{i+1}(perm);
end
subplot(121); hold on
histogram(mean(gScoresMat,2),'FaceColor','w','EdgeColor','k')
plot(mean(mean(gScoresMat,2))*ones(2,1),[0,max(get(gca,'YLim'))])
xlabel('mean score across nulls (should average to zero for no-bias)')
[~,ix] = sort(abs(mean(gScoresMat,2)),'descend');
subplot(122)
histogram(gScoresMat(ix(1),:),'FaceColor','w','EdgeColor','k');
title(sprintf('gene scores for %s',geneInfo.acronym{ix(1)}))
xlabel(sprintf('scores across %u nulls',numNulls))

%-------------------------------------------------------------------------------
% Produce some summary plots:
%-------------------------------------------------------------------------------
titleText = sprintf('%s-%s',whatEdgeMeasure,randomizeHow);



%-------------------------------------------------------------------------------
% Look at distribution for some top ones
%-------------------------------------------------------------------------------
f = figure('color','w');
for i = 1:min(15,numGOCategories)
    subplot(5,3,i); hold on
    histogram(categoryScores(ix_GO(i),nullInd),'edgeColor','k','FaceColor','w');
    plot(categoryScores(ix_GO(i),1)*ones(2,1),[0,max(get(gca,'ylim'))],'-r')
    % plot(whatStat(ix_GO(i))*ones(2,1),[0,max(get(gca,'ylim'))],'-r')
    title(sprintf('%s (%u; p_{corr}=%.2g)\n',GOTable.GOName{ix_GO(i)},...
                        sizeGOCategories(ix_GO(i)),pValsZ_corr(ix_GO(i))));
    % ,pValsZ(ix(i))
end
