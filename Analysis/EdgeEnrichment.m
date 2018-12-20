function [GOTable,gScore] = EdgeEnrichment(edgeData,nullEdgeData,whatNull,params)
% Computes correlation between a pairwise measure and gene expression outer
% product
%-------------------------------------------------------------------------------

% Processing parameters:
if nargin < 1
    whatEdgeMeasure = 'distance';
end
if nargin < 2
    nullEdgeData = [];
end
if nargin < 3
    whatNull = 'randomGene'; % 'topology', 'uniformTopology', 'permutedGeneDep', 'shuffleEdgeVals'
end
if nargin < 4
    params = GiveMeDefaultParams('mouse');
    % params.gcc.pValOrStat = 'stat'; % 'pval','stat'
    % params.gcc.thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
end

onlyOnEdges = params.gcc.onlyConnections; % only look where there are structural connections
corrType = params.gcc.whatCorr; % 'Spearman', 'Pearson'
correctDistance = params.gcc.regressDistance; % whether to regress distance
absType = params.gcc.absType;

% %-------------------------------------------------------------------------------
% % Checks:
% if strcmp(whatEdgeMeasure,'connected')
%     corrType = 'ttest';
% end
%
% %===============================================================================
% % Get data:
% [A_wei,regionAcronyms,A_p] = GiveMeAdj(params.c.connectomeSource,params.c.pThreshold,false,...
%                                     params.c.whatWeightMeasure,params.c.whatHemispheres,'all');


% % Make sure structInfo match with regionAcronyms:
% [~,ia,ib] = intersect(structInfo.acronym,regionAcronyms,'stable');
% if length(ia)~=height(structInfo) % this should be the limiting list (all should appear in regionAcronyms)
%     error('Error matching ROIs by acronym');
% end
% regionAcronyms = regionAcronyms(ib);
% A_wei = A_wei(ib,ib);
% if ~isempty(A_p)
%     A_p = A_p(ib,ib);
% end
% fprintf(1,'Rearranged information for %u structures\n',length(ib));
%
% % Filter structures:
% [A_wei,geneData,structInfo,keepInd] = filterStructures(params.c.structFilter,...
%                                                 structInfo,A_wei,geneData);
% A_bin = (A_wei~=0);
%
% if ~isempty(A_p)
%     A_p = A_p(keepInd,keepInd);
% end
%
% %-------------------------------------------------------------------------------
% % Compute the edge measure:
% %-------------------------------------------------------------------------------
% if strcmp(whatNull,'randomGene')
%     edgeData = GiveMeEdgeMeasure(whatEdgeMeasure,A_bin,A_wei,onlyOnEdges,whatSpecies,A_p);
% else
%     % Compute edge measures (+nulls):
%     edgeMeasures = GiveMeNullEdgeMeasures(whatNull,whatEdgeMeasure,A_bin,A_wei,numNulls,onlyOnEdges,structInfo);
% end

%-------------------------------------------------------------------------------
% Load in gene expression data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);

%-------------------------------------------------------------------------------
% Check that regions in the gene expression data match the edge data provided:
if size(geneData,1)~=length(edgeData)
    error('You loser, they''re not even the same size');
end

%-------------------------------------------------------------------------------
% Regress distance?:
if correctDistance
    distanceRegressor = GiveMeDistanceMatrix(params.humanOrMouse);
    fprintf(1,'Regressing inter-areal distances for %s\n',params.humanOrMouse);
else
    distanceRegressor = []; % just compute normal correlations
    fprintf(1,'No distance regressor used\n');
end

%-------------------------------------------------------------------------------
% Compute the scores -> GO enrichment
switch whatNull
case 'randomGene'
    % ---Conventional enrichment analysis under random-gene nulls---
    % (computational benefit: each gene scored independently allows
    % cheap computation of large number of null samples)
    fprintf(1,'Random gene nulls\n');

    % Compute gene scores:
    [geneScores,geneEntrezIDs] = GiveMeGCC(edgeData,geneData,geneInfo.entrez_id,params,distanceRegressor);

    % In-house category-size based random-gene permutation testing enrichment:
    GOTable = SingleEnrichment(geneScores,geneEntrezIDs,params.e);

    % Any significant?:
    numSig = sum(GOTable.pValCorr < params.e.sigThresh);
    fprintf(1,'%u significant categories at p_corr < %.1g\n',numSig,params.e.sigThresh);
    display(GOTable(1:numSig,:));

otherwise
    % ---Custom null---
    numNulls = params.gcc.numNulls;
    fprintf(1,'%u nulls\n',numNulls);

    % (include only annotations for genes with entrez IDs that are in our dataset)
    fprintf(1,'---Interested in %s---\n',params.e.processFilter);
    GOTable = GiveMeGOData(params,geneInfo.entrez_id);
    numGOCategories = height(GOTable);

    categoryScores = nan(numGOCategories,numNulls+1);
    geneScores = cell(numNulls+1,1);
    entrezIDsKept = cell(numNulls+1,1); % (sometimes entrez IDs change -- e.g., when matching to distance results)
    timer = tic;
    for i = 1:numNulls+1
        if i==1
            fprintf(1,'Read edge data!\n\n');
            theEdgeData = edgeData;
        else
            fprintf(1,'Null edge data: %u/%u\n\n',i-1,numNulls);
            theEdgeData = nullEdgeData{i-1};
        end

        % Compute gene scores:
        [geneScores{i},entrezIDsKept{i}] = GiveMeGCC(theEdgeData,geneData,geneInfo.entrez_id,params,distanceRegressor);

        % Record mean scores for each category:
        for j = 1:numGOCategories
            matchMe = ismember(entrezIDsKept{i},GOTable.annotations{j});
            if sum(matchMe) > 0
                categoryScores(j,i) = nanmean(geneScores{i}(matchMe));
            end
        end

        % Give user feedback (not so useful for parfor... :-/)
        % fprintf(1,'\n\n----%u/%u (%s remaining)\n\n',i,numNulls+1,...
        %                         BF_thetime((numNulls+1-i)*(toc(timer)/i)));
    end

    %-------------------------------------------------------------------------------
    % Compute p-values (+ corrected) and annotate GOTable
    GOTable = EstimatePVals(categoryScores,params.gcc.whatTail,GOTable);

    %-------------------------------------------------------------------------------
    % List categories with lowest p-values, or highest mean across nulls, etc.
    ix_GO = ListCategories(geneInfo,GOTable);
    GOTable = GOTable(ix_GO,:);
    categoryScores = categoryScores(ix_GO);

    %-------------------------------------------------------------------------------
    % Check whether the mean null score for each gene is zero:
    f = figure('color','w');
    allKeptEntrez = unique(vertcat(entrezIDsKept{:}));
    gScoresMat = nan(length(allKeptEntrez),numNulls);
    for i = 1:numNulls
        [~,~,perm] = intersect(allKeptEntrez,entrezIDsKept{i+1},'stable');
        gScoresMat(:,i) = geneScores{i+1}(perm);
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
    titleText = sprintf('%s-%s',whatEdgeMeasure,whatNull);
    NullSummaryPlots(GOTable,categoryScores,titleText);

    %-------------------------------------------------------------------------------
    % Look at distribution for some top ones:
    SpecificNullPlots(GOTable,categoryScores);
end

end
