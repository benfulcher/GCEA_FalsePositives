function ComputeRandomGeneNull(whatSpecies)
% Investigate how the distribution of correlations vary across all genes
%-------------------------------------------------------------------------------

whatSpecies = 'mouse';
%-------------------------------------------------------------------------------
% Bits of parameters:
params = GiveMeDefaultParams(whatSpecies);
numNullMaps = 250;
whatNullType = 'randomMap';
whatCorr = 'Spearman';

%-------------------------------------------------------------------------------
% Get real data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
numGenes = height(geneInfo);
numAreas = height(structInfo);
geneDataRand = rand(numAreas,numGenes);

%-------------------------------------------------------------------------------
% Under random phenotypes:
nullMaps = rand(numAreas,numNullMaps);
nullScores = cell(numNullMaps,1);
nullScoresRand = cell(numNullMaps,1);
for i = 1:numNullMaps
    fprintf(1,'Null map %u/%u\n',i,numNullMaps);
    myNullMap = nullMaps(:,i);

    % Correlation of null map to every real gene:
    nullScoresHere = nan(numGenes,1);
    parfor j = 1:numGenes
        nullScoresHere(j) = corr(myNullMap,geneData(:,j),'type',whatCorr,'rows','pairwise');
    end
    nullScores{i} = nullScoresHere;

    % Correlation of null map to every random gene:
    nullScoresRandomizedGeneData = nan(numGenes,1);
    parfor j = 1:numGenes
        nullScoresRandomizedGeneData(j) = corr(myNullMap,geneDataRand(:,j),'type',whatCorr,'rows','pairwise');
    end
    nullScoresRand{i} = nullScoresRandomizedGeneData;
end

nullScoreMeans = {cellfun(@mean,nullScores),cellfun(@mean,nullScoresRand)};
nullScoreVar = {cellfun(@var,nullScores),cellfun(@var,nullScoresRand)};

%-------------------------------------------------------------------------------
f = figure('color','w');
extraParams = struct();
extraParams.theColors = BF_getcmap('set1',3,1);
extraParams.customSpot = '';
extraParams.offsetRange = 0.7;
BF_JitteredParallelScatter(nullScoreMeans,true,true,false,extraParams);

f = figure('color','w');
extraParams = struct();
extraParams.theColors = BF_getcmap('set1',3,1);
extraParams.customSpot = '';
extraParams.offsetRange = 0.7;
BF_JitteredParallelScatter(nullScoreVar,true,true,false,extraParams);

f = figure('color','w');
hold('on')
plot(nullScoreMeans{1},nullScoreVar{1},'.k')
plot(nullScoreMeans{2},nullScoreVar{2},'.r')
xlabel('mean')
ylabel('variance')
title('Properties of the distribution of correlation coefficients across all genes to different null maps')

end
