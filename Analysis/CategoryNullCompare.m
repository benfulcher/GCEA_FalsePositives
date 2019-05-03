%-------------------------------------------------------------------------------
% CategoryNullCompare
% Compares different nulls for correlations with some underlying phenotype
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Settings:
%-------------------------------------------------------------------------------
whatGOID = 6099;
whatShuffle = 'twoIsocortex';
whatSpecies = 'mouse';
whatCorr = 'Spearman'; % 'Pearson', 'Spearman'
numNulls = 100;
whatPhenotype = 'degree';

%-------------------------------------------------------------------------------
% Load in gene-expression data for this GO category:
params = GiveMeDefaultParams(whatSpecies);
[geneData,geneInfo,structInfo,categoryInfo] = GiveMeGOCategory(whatGOID,params);
numGenesGO = height(geneInfo);

%-------------------------------------------------------------------------------
% Load the phenotype map:
switch whatPhenotype
case 'degree'
    % Load adjacency matrix:
    [A_bin,regionAcronyms,adjPVals] = GiveMeAdj(params.c.connectomeSource,params.c.pThreshold,true,...
                                        params.c.whatWeightMeasure,params.c.whatHemispheres,...
                                        params.c.structFilter);
    thePhenotype = sum(A_bin,1)' + sum(A_bin,2);
end

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
    % Coherent spatial shuffling:
    %-------------------------------------------------------------------------------
    gScore = zeros(numGenesGO,2);
    for i = 1:numGenesGO
        gScore(i,1) = corr(thePhenotype,geneData(permVectorAll,i),'type',whatCorr,'rows','pairwise');
        gScore(i,2) = corr(thePhenotype,geneData(permVectorCustom,i),'type',whatCorr,'rows','pairwise');
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
            gScore(i) = corr(thePhenotype,geneData(:,rp(i)),'type',whatCorr,'rows','pairwise');
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
title(sprintf('%s\t%s (%u)',categoryInfo.GOIDlabel{1},...
                    categoryInfo.GOName{1},categoryInfo.size))
xlabel(sprintf('%s correlation',whatCorr))
legend([h1,h2,h3],{'random-area-all',sprintf('random-area-%s',whatShuffle),'random-gene'})
