%-------------------------------------------------------------------------------
% CategoryNullCompare
% Does nulls for correlations with degree, in different ways
%-------------------------------------------------------------------------------

% Settings:
whatGOID = 6099;
whatCorr = 'Spearman'; % 'Pearson', 'Spearman'
numNulls = 2000;
whatShuffle = 'twoIsocortex';
whatSpecies = 'mouse';

%-------------------------------------------------------------------------------
% Load defaults:
params = GiveMeDefaultParams(whatSpecies);
warning('Over-writing gene normalization settings -> none')
params.g.normalizationGene = 'none';
params.g.normalizationRegion = 'none';

% Get gene data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);

% Load in GO annotations:
GOTable = GetFilteredGOData(params.e.whatSource,params.e.processFilter,params.e.sizeFilter,...
                                    geneInfo.entrez_id);
numGOCategories = height(GOTable);
whatCategory = find(GOTable.GOID==whatGOID);
fprintf(1,'Looking in at %s:%s (%u)\n',GOTable.GOIDlabel{whatCategory},...
                    GOTable.GOName{whatCategory},GOTable.size(whatCategory));

% Load adjacency matrix:
[A_bin,regionAcronyms,adjPVals] = GiveMeAdj(params.c.connectomeSource,params.c.pThreshold,true,...
                                    params.c.whatWeightMeasure,params.c.whatHemispheres,...
                                    params.c.structFilter);
k = sum(A_bin,1)' + sum(A_bin,2);

theGenesEntrez = GOTable.annotations{whatCategory};
[entrezMatched,ia,ib] = intersect(theGenesEntrez,geneInfo.entrez_id);
fprintf(1,'%u/%u genes from this GO category match\n',length(entrezMatched),length(theGenesEntrez));
numGenesGO = length(entrezMatched);
matchMe = find(ismember(geneInfo.entrez_id,entrezMatched));

categoryScores = nan(numNulls+1,3);
for n = 1:numNulls+1

    if n == 1
        fprintf(1,'True ordering!');
        permVectorAll = 1:size(geneData,1);
        permVectorCustom = 1:size(geneData,1);
    else
        fprintf(1,'Null %u/%u\n',n,numNulls+1);
        permVectorAll = randperm(size(geneData,1));
        permVectorCustom = AnatomyShuffle(structInfo.divisionLabel,whatShuffle);
    end

    %-------------------------------------------------------------------------------
    % Shuffling node properties:
    %-------------------------------------------------------------------------------
    gScore = zeros(numGenesGO,2);
    for i = 1:numGenesGO
        gScore(i,1) = corr(k,geneData(permVectorAll,matchMe(i)),'type',whatCorr,'rows','pairwise');
        gScore(i,2) = corr(k,geneData(permVectorCustom,matchMe(i)),'type',whatCorr,'rows','pairwise');
    end

    % Record mean scores for each category:
    categoryScores(n,1) = nanmean(gScore(:,1));
    categoryScores(n,2) = nanmean(gScore(:,2));

    %-------------------------------------------------------------------------------
    % Taking genes at random:
    %-------------------------------------------------------------------------------
    if n == 1
        categoryScores(1,3) = categoryScores(1,1);
    else
        rp = randperm(size(geneData,2),numGenesGO);
        gScore = zeros(numGenesGO,1);
        for i = 1:numGenesGO
            gScore(i) = corr(k,geneData(:,rp(i)),'type',whatCorr,'rows','pairwise');
        end
        categoryScores(n,3) = nanmean(gScore);
    end

    %-------------------------------------------------------------------------------
    % Both randomized!:
    %-------------------------------------------------------------------------------
    % rp = randperm(size(geneData,2),length(theGenesEntrez));
    % gScore = zeros(length(theGenesEntrez),1);
    % for i = 1:length(theGenesEntrez)
    %     gScore(i) = corr(k,geneData(permVector,rp(i)),'type',whatCorr,'rows','pairwise');
    % end
    % categoryScores(n,3) = nanmean(gScore);
end

%-------------------------------------------------------------------------------
% Estimate p-values
%-------------------------------------------------------------------------------
nullInd = 2:numNulls+1;
pValZ = arrayfun(@(x)1-normcdf(categoryScores(1,x),mean(categoryScores(nullInd,x)),...
                    std(categoryScores(nullInd,x))),1:3);
labels = {'shuffle-all',whatShuffle,'random-gene'};
for i = 1:3
    fprintf(1,'(%s): p_Z = %g\n',labels{i},pValZ(i));
end

%-------------------------------------------------------------------------------
% Plot null distributions:
%-------------------------------------------------------------------------------
f = figure('color','w'); ax = gca; hold on
% First the three nulls as histograms:
h1 = histogram(categoryScores(2:end,1)); h1.FaceColor = [0.2157    0.4941    0.7216];
h2 = histogram(categoryScores(2:end,2)); h2.FaceColor = [0.5961    0.3059    0.6392];
h3 = histogram(categoryScores(2:end,3)); h3.FaceColor = [0.3020    0.6863    0.2902];
plot(ones(2,1)*mean(categoryScores(2:end,1)),ax.YLim,'color',h1.FaceColor);
plot(ones(2,1)*mean(categoryScores(2:end,2)),ax.YLim,'color',h2.FaceColor);
plot(ones(2,1)*mean(categoryScores(2:end,3)),ax.YLim,'color',h3.FaceColor);
plot(ones(2,1)*categoryScores(1,1),ax.YLim,'r')
title(sprintf('%s\t%s (%u)',GOTable.GOIDlabel{whatCategory},...
                    GOTable.GOName{whatCategory},GOTable.size(whatCategory)))
xlabel(sprintf('%s correlation',whatCorr))
legend([h1,h2,h3],{'random-area-all',sprintf('random-area-%s',whatShuffle),'random-gene'})

%-------------------------------------------------------------------------------
f = figure('color','w'); ax = gca;
corrMat = corr(geneData(:,matchMe),'type',whatCorr,'rows','pairwise');
ord = BF_ClusterReorder(corrMat,'euclidean','average');
imagesc(corrMat(ord,ord));
ax.YTick = 1:length(matchMe);
acronymsMatch = geneInfo.acronym(matchMe);
ax.YTickLabel = acronymsMatch(ord);
caxis([-1,1])
colormap([flipud(BF_getcmap('blues',9));1,1,1;BF_getcmap('reds',9)])
title(sprintf('%s %s (%u genes)[%s]',GOTable.GOIDlabel{whatCategory},...
                    GOTable.GOName{whatCategory},GOTable.size(whatCategory),...
                    whatCorr))
axis('square')
cB = colorbar;
