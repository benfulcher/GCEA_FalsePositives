function [GOTable,gScore,geneEntrezIDs] = EdgeEnrichment(whatEdgeMeasure,...
                                            onlyOnEdges,correctDistance,absType,corrType)
% Computes correlation between a pairwise measure and gene expression outer
% product
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Fixed processing parameters:
%-------------------------------------------------------------------------------
if nargin < 1
    whatEdgeMeasure = 'distance';
end
if nargin < 2
    onlyOnEdges = false;
end
if nargin < 3
    correctDistance = false; % false,true;
end
if nargin < 4
    absType = 'pos'; % 'pos','neg','abs' -> e.g., pos -> coexpression contribution increases with the statistic
end
if nargin < 5
    corrType = 'Spearman'; % {'Spearman','Pearson'};
end

%===============================================================================
% Connectome processing
cParam = GiveMeDefaultParams('conn');

% Gene processing
gParam = GiveMeDefaultParams('gene');
gParam.subsetOfGenes = []; % only look at the first X genes. Set to empty for all
                         % genes & save a .mat file of results

% Null model:
whatNull = 'randomGene'; % 'topology', 'uniformTopology', 'permutedGeneDep', 'shuffleEdgeVals'
numNulls = 100;

% Correlations:
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
pValOrStat = 'stat'; % 'pval','stat'

% Enrichment parameters:
eParam = GiveMeDefaultParams('enrichment');

%-------------------------------------------------------------------------------
% Checks:
if strcmp(whatEdgeMeasure,'connected')
    corrType = 'ttest';
end

%===============================================================================
% Get data:
[A_wei,regionAcronyms,A_p] = GiveMeAdj(cParam.connectomeSource,cParam.pThreshold,false,...
                                    cParam.whatWeightMeasure,cParam.whatHemispheres,'all');
A_bin = (A_wei~=0);
[geneData,geneInfo,structInfo] = LoadMeG({gParam.normalizationGene,gParam.normalizationRegion},gParam.energyOrDensity);
[A_bin,geneData,structInfo,keepInd] = filterStructures(cParam.structFilter,structInfo,A_bin,geneData);
A_wei = A_wei(keepInd,keepInd);
A_p = A_p(keepInd,keepInd);

%-------------------------------------------------------------------------------
% Compute the edge measure:
%-------------------------------------------------------------------------------
if strcmp(whatNull,'randomGene')
    edgeData = GiveMeEdgeMeasure(whatEdgeMeasure,A_bin,A_wei,onlyOnEdges,A_p);
else
    % Compute edge measures (+nulls):
    edgeMeasures = GiveMeNullEdgeMeasures(whatNull,whatEdgeMeasure,A_bin,A_wei,numNulls,onlyOnEdges);
end

%-------------------------------------------------------------------------------
% Regress distance?:
%-------------------------------------------------------------------------------
if correctDistance
    C = load('Mouse_Connectivity_Data.mat','Dist_Matrix');
    distanceRegressor = C.Dist_Matrix{1,1};
    fprintf(1,'Regressing ipsilateral distances\n');
else
    fprintf(1,'No distance regressor used\n');
    distanceRegressor = []; % just compute normal correlations
end

%-------------------------------------------------------------------------------
% Compute the scores -> GO enrichment
%-------------------------------------------------------------------------------
switch whatNull
case 'randomGene'
    % Conventional gene enrichment analysis
    %-------------------------------------------------------------------------------
    % Compute gene scores:
    %-------------------------------------------------------------------------------
    [gScore,geneEntrezIDs] = GiveMeGCC(edgeData,geneData,geneInfo.entrez_id,corrType,...
                            distanceRegressor,absType,thresholdGoodGene,pValOrStat);

    %-------------------------------------------------------------------------------
    % Enrichment using our in-house random-gene null method:
    %-------------------------------------------------------------------------------
    [GOTable,geneEntrezAnnotations] = SingleEnrichment(gScore,geneEntrezIDs,...
                        eParam.processFilter,eParam.sizeFilter,eParam.numIterations);

    %-------------------------------------------------------------------------------
    % ANALYSIS:
    %-------------------------------------------------------------------------------
    numSig = sum(GOTable.pVal_corr < eParam.enrichmentSigThresh);
    fprintf(1,'%u significant categories at p_corr < %.1g\n',numSig,eParam.enrichmentSigThresh);
    display(GOTable(1:numSig,:));

otherwise
    %-------------------------------------------------------------------------------
    % Get GO data
    % (include only annotations for genes with entrez IDs that are in our dataset)
    [GOTable,geneEntrezAnnotations] = GetFilteredGOData(processFilter,sizeFilter,geneInfo.entrez_id);
    sizeGOCategories = cellfun(@length,geneEntrezAnnotations);
    numGOCategories = height(GOTable);

    categoryScores = nan(numGOCategories,numNulls+1);
    gScores = cell(numNulls+1,1);
    entrezIDsKept = cell(numNulls+1,1);
    fprintf(1,'---Interested in Biological Processes with FDR p < %g\n',eParam.enrichmentSigThresh);
    timer = tic;
    parfor i = 1:numNulls+1
        fprintf(1,'%u/%u\n\n',i,numNulls+1);

        % Compute gene scores:
        % (sometimes entrez IDs change -- e.g., when matching to distance results)
        theEdgeData = edgeMeasures{i};

        % Compute the score:
        [gScores{i},entrezIDsKept{i}] = GiveMeGCC(theEdgeData,theGeneData,...
                            geneInfo.entrez_id,corrType,distanceRegressor,absType,...
                            thresholdGoodGene,pValOrStat);

        % Record mean scores for each category:
        for j = 1:numGOCategories
            matchMe = ismember(entrezIDsKept{i},geneEntrezAnnotations{j});
            if sum(matchMe) <= 1
                continue
            end
            categoryScores(j,i) = nanmean(gScores{i}(matchMe));
        end

        % Give user feedback
        % fprintf(1,'\n\n----%u/%u (%s remaining)\n\n',i,numNulls+1,...
        %                         BF_thetime((numNulls+1-i)*(toc(timer)/i)));
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
    NullSummaryPlots(pValsZ,pValsZ_corr,categoryScores,meanNull,stdNull,sizeGOCategories,titleText);

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
end



end
