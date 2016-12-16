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
whatEdgeMeasure = 'bin_edgeBet';
onlyOnEdges = true; % whether to put values only on existing edges
                    % (rather than all node pairs for some measures)

% Randomization
randomizeHow = 'permutedGeneDep'; % 'uniformTopology', 'permutedGeneDep'
numNulls = 250;

% Gene processing
energyOrDensity = 'energy'; % what gene expression data to use
normalizationGene = 'none'; % 'none', 'mixedSigmoid'
normalizationRegion = 'none'; % 'none', 'zscore'

% GO settings
processFilter = 'biological_process';
sizeFilter = [5,200];

% Correlations:
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
corrType = 'Spearman'; % {'Spearman','Pearson'};
pValOrStat = 'stat'; % 'pval','stat'
absType = 'pos'; % 'pos','neg','abs' -> e.g., pos -> coexpression contribution increases with the statistic
correctDistance = false; % false,true;


%-------------------------------------------------------------------------------
% Define a set of edge measures to compare:
[A_bin,regionStruct,adjPVals] = GiveMeAdj(connectomeSource,pThreshold,true,whatHemispheres,justCortex);
A_wei = GiveMeAdj(connectomeSource,pThreshold,false,whatHemispheres,justCortex);

%-------------------------------------------------------------------------------
% Compute randomized
switch randomizeHow
case 'uniformTopology'
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

        % weightVector = A_wei(A_wei > 0);
        % A_wei_rand = A_rand;
        % A_wei_rand(A_wei_rand>0) = weightVector;

        switch whatEdgeMeasure
        case 'ktotktot'
            ktot = sum(A_rand,1)' + sum(A_rand,2);
            product = ktot*ktot';
            product(A_rand == 0) = 0;
            edgeMeasures{i} = product;
        case 'bin_edgeBet'
            edgeMeasures{i} = edge_betweenness_bin(A_rand);
        end
    end
case 'permutedGeneDep'
    edgeMeasures = edge_betweenness_bin(A_bin);
end

%-------------------------------------------------------------------------------
% Retrieve and match gene data
[geneData,geneInfo,structInfo] = LoadMeG({normalizationGene,normalizationRegion},energyOrDensity);
% Check gene data matches connectome data
if ~all([regionStruct.id]'==structInfo.id)
    % Take subset
    [~,ia,ib] = intersect([regionStruct.id]',structInfo.id,'stable');
    geneData = geneData(ib,:);
    structInfo = structInfo(ib,:);
    fprintf(1,'Gene data matched to subset of %u Allen regions\n',length(ib));
end
numGenes = size(geneData,1);
% Don't regress distance:
distanceRegressor = [];

%-------------------------------------------------------------------------------
% Get GO data
[GOTable,geneEntrezAnnotations] = GetFilteredGOData(processFilter,sizeFilter);
sizeGOCategories = cellfun(@length,geneEntrezAnnotations);
numGOCategories = height(GOTable);

%-------------------------------------------------------------------------------
% Get scores for genes:
%-------------------------------------------------------------------------------
categoryScores = nan(numGOCategories,numNulls+1);
enrichmentSigThresh = 0.05;
fprintf(1,'---Interested in Biological Processes with FDR p < %g\n',enrichmentSigThresh);
timer = tic;
for i = 1:numNulls+1
    fprintf(1,'%u/%u\n\n',i,numNulls+1);

    % Compute gene scores:
    % (sometimes entrez IDs change -- e.g., when matching to distance results)
    switch randomizeHow
    case 'uniformTopology'
        % each null has a different set of edge measures:
        theEdgeData = edgeMeasures{i};
        theGeneData = geneData;
    case 'permutedGeneDep'
        % each null should use the same edge measure but different permutations of the gene data
        theEdgeData = edgeMeasures;
        rp = randperm(numGenes);
        theGeneData = geneData(rp,:);
    end
    % Compute the score:
    [gScore,geneEntrezIDs] = GiveMeGCC(theEdgeData,theGeneData,...
                            geneInfo.entrez_id,corrType,distanceRegressor,absType,...
                            thresholdGoodGene,pValOrStat);

    % Record mean scores for each category:
    for j = 1:numGOCategories
        matchMe = ismember(geneEntrezIDs,geneEntrezAnnotations{j});
        if sum(matchMe) <= 1
            continue
        end
        categoryScores(j,i) = nanmean(gScore(matchMe));
    end

    % Give user feedback
    fprintf(1,'\n\n----%u/%u (%s remaining)\n\n',i,numNulls+1,...
                            BF_thetime((numNulls+1-i)*(toc(timer)/i)));
end

%-------------------------------------------------------------------------------
% Save
%-------------------------------------------------------------------------------
fileName = sprintf('%s-%s-%s-%unulls.mat',whatEdgeMeasure,randomizeHow,processFilter,numNulls);
save(fullfile('DataOutputs',fileName));

%-------------------------------------------------------------------------------
% List categories with the highest mean nulls
%-------------------------------------------------------------------------------
numTop = 100;
nullInd = 2:numNulls+1;
meanNull = nanmean(categoryScores(:,nullInd),2); % mean score of genes in each category
stdNull = nanstd(categoryScores(:,nullInd),[],2); % std of genes in each category
% looking at right tail:
pValsPerm = arrayfun(@(x)mean(categoryScores(x,nullInd)>=categoryScores(x,1)),1:numGOCategories);
pValsZ = arrayfun(@(x)1-normcdf(categoryScores(x,1),mean(categoryScores(x,nullInd)),std(categoryScores(x,nullInd))),1:numGOCategories);
pValsZ_corr = mafdr(pValsZ,'BHFDR','true');
whatStat = pValsZ_corr;
[~,ix] = sort(whatStat,'ascend');
fprintf(1,'%u nans removed\n',sum(isnan(whatStat)));
ix(isnan(whatStat(ix))) = [];
for i = 1:numTop
    fprintf(1,'%u (%u genes): %s (%.2g)\n',i,sizeGOCategories(ix(i)),...
                        GOTable.GOName{ix(i)},whatStat(ix(i)));
end

% Plot distribution of p-values:
f = figure('color','w'); hold on
histogram(pValsZ)
histogram(pValsZ_corr)
xlabel('p-values')
legend({'raw','corrected'})
ylabel('frequency')

% Plot distribution of mean nulls:
f = figure('color','w');
histogram(meanNull)
xlabel('mean corr statistic across nulls')
ylabel('frequency')

% Check dependence on GO category size
f = figure('color','w');
plot(sizeGOCategories,whatStat,'.k')
xlabel('GO category size')
ylabel('corrected p-value')

% Relationship between null variance and corrected p-value
f = figure('color','w'); hold on
plot(meanNull,categoryScores(:,1),'.k')
plot([min(meanNull),max(meanNull)],[min(meanNull),max(meanNull)],'r')
xlabel('mean of null distribution')
ylabel('real scores')

% Std
f = figure('color','w'); hold on
plot(stdNull,categoryScores(:,1),'.k')
xlabel('std of null distribution')
ylabel('real scores')

%-------------------------------------------------------------------------------
% Look at distribution for some top ones
%-------------------------------------------------------------------------------
f = figure('color','w');
for i = 1:10
    subplot(5,2,i); hold on
    histogram(categoryScores(ix(i),nullInd));
    plot(categoryScores(ix(i),1)*ones(2,1),[0,max(get(gca,'ylim'))],'-r')
    % plot(whatStat(ix(i))*ones(2,1),[0,max(get(gca,'ylim'))],'-r')
    title(GOTable.GOName{ix(i)})
    title(sprintf('%s (%u genes; p = %.2g)\n',GOTable.GOName{ix(i)},...
                        sizeGOCategories(ix(i)),pValsZ_corr(ix(i))));
end
