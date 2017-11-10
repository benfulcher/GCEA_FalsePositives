function [GOTable,gScore] = EdgeEnrichment(whatEdgeMeasure,onlyOnEdges,...
                correctDistance,absType,corrType,whatNull,numNulls,whatSpecies,params)
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
if nargin < 6
    whatNull = 'randomGene'; % 'topology', 'uniformTopology', 'permutedGeneDep', 'shuffleEdgeVals'
end
if nargin < 7
    numNulls = 100; % not used for 'randomGene'
end
if nargin < 8
    whatSpecies = 'mouse';
    fprintf(1,'Mouse by default\n');
end
if nargin < 9
    % Just use global defaults if none provided
    params = GiveMeDefaultParams(whatSpecies);
end

%-------------------------------------------------------------------------------
% Settings for computing correlations:
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
pValOrStat = 'stat'; % 'pval','stat'

%-------------------------------------------------------------------------------
% Checks:
if strcmp(whatEdgeMeasure,'connected')
    corrType = 'ttest';
end

%===============================================================================
% Get data:
[A_wei,regionAcronyms,A_p] = GiveMeAdj(params.c.connectomeSource,params.c.pThreshold,false,...
                                    params.c.whatWeightMeasure,params.c.whatHemispheres,'all');
[geneData,geneInfo,structInfo] = LoadMeG(params.g);

% Make sure structInfo match with regionAcronyms:
[~,ia,ib] = intersect(structInfo.acronym,regionAcronyms,'stable');
if length(ia)~=height(structInfo) % this should be the limiting list (all should appear in regionAcronyms)
    error('Error matching ROIs by acronym');
end
regionAcronyms = regionAcronyms(ib);
A_wei = A_wei(ib,ib);
if ~isempty(A_p)
    A_p = A_p(ib,ib);
end
fprintf(1,'Rearranged information for %u structures\n',length(ib));

% Filter structures:
[A_wei,geneData,structInfo,keepInd] = filterStructures(params.c.structFilter,...
                                                structInfo,A_wei,geneData);
A_bin = (A_wei~=0);

if ~isempty(A_p)
    A_p = A_p(keepInd,keepInd);
end

%-------------------------------------------------------------------------------
% Compute the edge measure:
%-------------------------------------------------------------------------------
if strcmp(whatNull,'randomGene')
    edgeData = GiveMeEdgeMeasure(whatEdgeMeasure,A_bin,A_wei,onlyOnEdges,whatSpecies,A_p);
else
    % Compute edge measures (+nulls):
    edgeMeasures = GiveMeNullEdgeMeasures(whatNull,whatEdgeMeasure,A_bin,A_wei,numNulls,onlyOnEdges,structInfo);
end

%-------------------------------------------------------------------------------
% Regress distance?:
%-------------------------------------------------------------------------------
if correctDistance
    distanceRegressor = GiveMeDistanceMatrix(whatSpecies);
    fprintf(1,'Regressing distances\n');
else
    distanceRegressor = []; % just compute normal correlations
    fprintf(1,'No distance regressor used\n');
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
    GOTable = SingleEnrichment(gScore,geneEntrezIDs,...
                        params.e.whatSource,params.e.processFilter,params.e.sizeFilter,params.e.numIterations);

    %-------------------------------------------------------------------------------
    % ANALYSIS:
    %-------------------------------------------------------------------------------
    numSig = sum(GOTable.pValCorr < params.e.enrichmentSigThresh);
    fprintf(1,'%u significant categories at p_corr < %.1g\n',numSig,params.e.enrichmentSigThresh);
    display(GOTable(1:numSig,:));

otherwise
    %-------------------------------------------------------------------------------
    % Get GO data
    % (include only annotations for genes with entrez IDs that are in our dataset)
    GOTable = GetFilteredGOData(params.e.whatSource,Param.processFilter,...
                                    params.e.sizeFilter,geneInfo.entrez_id);
    numGOCategories = height(GOTable);

    categoryScores = nan(numGOCategories,numNulls+1);
    gScore = cell(numNulls+1,1);
    entrezIDsKept = cell(numNulls+1,1);
    fprintf(1,'---Interested in Biological Processes with FDR p < %g\n',...
                        params.e.enrichmentSigThresh);
    timer = tic;
    parfor i = 1:numNulls+1
        fprintf(1,'%u/%u\n\n',i,numNulls+1);

        % Compute gene scores:
        % (sometimes entrez IDs change -- e.g., when matching to distance results)
        [gScore{i},entrezIDsKept{i}] = GiveMeGCC(edgeMeasures{i},geneData,...
                            geneInfo.entrez_id,corrType,distanceRegressor,absType,...
                            thresholdGoodGene,pValOrStat);

        % Record mean scores for each category:
        for j = 1:numGOCategories
            matchMe = ismember(entrezIDsKept{i},GOTable.annotations{j});
            if sum(matchMe) <= 1
                continue
            end
            categoryScores(j,i) = nanmean(gScore{i}(matchMe));
        end

        % Give user feedback (not so useful for parfor... :-/)
        % fprintf(1,'\n\n----%u/%u (%s remaining)\n\n',i,numNulls+1,...
        %                         BF_thetime((numNulls+1-i)*(toc(timer)/i)));
    end

    %-------------------------------------------------------------------------------
    % Compute p-values (+ corrected) and annotate GOTable
    %-------------------------------------------------------------------------------
    fprintf(1,'GO categories with correlations to %s than %s nulls\n',whatEdgeMeasure,whatNull);
    whatTail = 'right'; % right-tailed p-values
    GOTable = EstimatePVals(categoryScores,whatTail,GOTable);

    %-------------------------------------------------------------------------------
    % List categories with lowest p-values, or highest mean across nulls, etc.
    %-------------------------------------------------------------------------------
    ix_GO = ListCategories(geneInfo,GOTable);
    GOTable = GOTable(ix_GO,:);

    %-------------------------------------------------------------------------------
    % Save to mat file:
    %-------------------------------------------------------------------------------
    if isempty(params.g.subsetOfGenes)
        fileName = sprintf('%s-%s-%s-G%s_R%s-%unulls.mat',whatEdgeMeasure,whatNull,...
                params.e.processFilter,params.g.normalizationGene,params.g.normalizationRegion,numNulls);
        save(fullfile('DataOutputs',fileName));
        fprintf(1,'Saved %s\n',fileName);
    end

    %-------------------------------------------------------------------------------
    % Check whether the mean null score for each gene is zero:
    f = figure('color','w');
    allKeptEntrez = unique(vertcat(entrezIDsKept{:}));
    gScoresMat = nan(length(allKeptEntrez),numNulls);
    for i = 1:numNulls
        [~,~,perm] = intersect(allKeptEntrez,entrezIDsKept{i+1},'stable');
        gScoresMat(:,i) = gScore{i+1}(perm);
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
    SpecificNullPlots(categoryScores,GOTable,ix_GO);
end

end
