%-------------------------------------------------------------------------------
% CategoryNullCompare
%-------------------------------------------------------------------------------

whatCategory = 17;
whatCorr = 'Pearson';

% [GOTable,geneEntrezAnnotations] = GetFilteredGOData('biological_process',[5,100],geneInfo.entrez_id);
sizeGOCategories = cellfun(@length,geneEntrezAnnotations);
numGOCategories = height(GOTable);
k = sum(A_bin,1)' + sum(A_bin,2);

theGenesEntrez = geneEntrezAnnotations{whatCategory};
matchMe = find(ismember(geneInfo.entrez_id,theGenesEntrez));

numNulls = 500;
categoryScores = nan(numNulls+1,3);
for n = 1:numNulls+1
    fprintf(1,'Null %u/%u\n',n,numNulls+1);

    if n == 1
        permVector = 1:size(geneData,1);
    else
        % permVector = randperm(size(geneData,1));
        permVector = AnatomyShuffle(structInfo.divisionLabel,'twoBroad');
    end

    %-------------------------------------------------------------------------------
    % Shuffling node properties:
    %-------------------------------------------------------------------------------
    gScore = zeros(length(theGenesEntrez),1);
    for i = 1:length(theGenesEntrez)
        gScore(i) = corr(k,geneData(permVector,matchMe(i)),'type',whatCorr,'rows','pairwise');
    end

    % Record mean scores for each category:
    categoryScores(n,1) = nanmean(gScore);

    %-------------------------------------------------------------------------------
    % Taking genes at random:
    %-------------------------------------------------------------------------------
    rp = randperm(size(geneData,2),length(theGenesEntrez));
    gScore = zeros(length(theGenesEntrez),1);
    for i = 1:length(theGenesEntrez)
        gScore(i) = corr(k,geneData(:,rp(i)),'type',whatCorr,'rows','pairwise');
    end
    categoryScores(n,2) = nanmean(gScore);

    %-------------------------------------------------------------------------------
    % Both randomized!:
    %-------------------------------------------------------------------------------
    rp = randperm(size(geneData,2),length(theGenesEntrez));
    gScore = zeros(length(theGenesEntrez),1);
    for i = 1:length(theGenesEntrez)
        gScore(i) = corr(k,geneData(permVector,rp(i)),'type',whatCorr,'rows','pairwise');
    end
    categoryScores(n,3) = nanmean(gScore);

end

%-------------------------------------------------------------------------------
f = figure('color','w'); ax = gca; hold on
plot(ones(2,1)*mean(categoryScores(2:end,1)),ax.YLim,'color',h1.FaceColor);
plot(ones(2,1)*mean(categoryScores(:,3)),ax.YLim,'color',h3.FaceColor);
h1 = histogram(categoryScores(2:end,1)); h1.FaceColor = [0.2157    0.4941    0.7216];
h3 = histogram(categoryScores(:,3)); h3.FaceColor = [0.5961    0.3059    0.6392];
h2 = histogram(categoryScores(:,2)); h2.FaceColor = [0.3020    0.6863    0.2902];
plot(ones(2,1)*mean(categoryScores(:,2)),ax.YLim,'color',h2.FaceColor);
plot(ones(2,1)*categoryScores(1,1),ax.YLim,'r')
title(sprintf('%s\t%s (%u)',GOTable.GOIDlabel{whatCategory},GOTable.GOName{whatCategory},GOTable.size(whatCategory)))
xlabel(sprintf('%s correlation',whatCorr))
legend([h1,h2,h3],{'random-degree','random-gene','random-both'})

%-------------------------------------------------------------------------------
f = figure('color','w'); ax = gca;
corrMat = corr(geneData(:,matchMe),'type','Pearson','rows','pairwise');
ord = BF_ClusterReorder(corrMat,'euclidean','average')
imagesc(corrMat(ord,ord));
ax.YTick = 1:length(matchMe);
acronymsMatch = geneInfo.acronym(matchMe);
ax.YTickLabel = acronymsMatch(ord);
caxis([-1,1])
colormap([flipud(BF_getcmap('blues',9));1,1,1;BF_getcmap('reds',9)])
title(sprintf('%s\t%s (%u)',GOTable.GOIDlabel{whatCategory},GOTable.GOName{whatCategory},GOTable.size(whatCategory)))
axis('square')
