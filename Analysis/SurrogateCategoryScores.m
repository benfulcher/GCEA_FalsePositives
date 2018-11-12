%-------------------------------------------------------------------------------
% SurrogateCategoryScores
%-------------------------------------------------------------------------------
% Aim is to get a distribution of mean category scores for surrogate spatial maps
%-------------------------------------------------------------------------------

numMaps = 500;
whatSpecies = 'mouse';

%-------------------------------------------------------------------------------
% Prepare for enrichment:
params = GiveMeDefaultParams(whatSpecies);
switch whatSpecies
case 'mouse'
    realAndFake = {'mouse','surrogate-mouse'};
case 'human'
    realAndFake = {'human','surrogate-human'};
end

%-------------------------------------------------------------------------------
% Load data:
params.g.humanOrMouse = realAndFake{1};
[geneDataReal,geneInfoReal,structInfoReal] = LoadMeG(params.g);
params.g.humanOrMouse = realAndFake{2};
[geneDataNull,geneInfoNull,structInfoNull] = LoadMeG(params.g);
GOTable = GetFilteredGOData(params.e.dataSource,params.e.processFilter,...
                                params.e.sizeFilter,geneInfoReal.entrez_id);
numGenes = height(geneInfoReal);
numGOCategories = height(GOTable);
if numMaps > height(geneInfoNull)
    error('There aren''t enough null maps to compare against...');
end

%-------------------------------------------------------------------------------
% Assign scores to categories of genes (repeat across an ensemble of random spatially correlated maps)
categoryScoresRaw = nan(numGOCategories,numMaps);
categoryScoresAbs = nan(numGOCategories,numMaps);
parfor n = 1:numMaps
    map_n = geneDataNull(:,n);

    % Score all genes for spatial correlation with this map:
    gScore = zeros(numGenes,1);
    for i = 1:numGenes
        gScore(i) = corr(map_n,geneDataReal(:,i),'type','Spearman');
    end

    % Record mean scores for each category:
    for j = 1:numGOCategories
        matchMe = ismember(geneInfoReal.entrez_id,GOTable.annotations{j});
        if any(matchMe)
            categoryScoresRaw(j,n) = nanmean(gScore(matchMe));
            categoryScoresAbs(j,n) = nanmean(abs(gScore(matchMe)));
        end
    end
end

%-------------------------------------------------------------------------------
% Find categories with highest scores
meanScore = nanmean(squeeze(categoryScoresAbs),2);
GOTable.meanAbsScore = meanScore;
[~,ix] = sort(meanScore,'descend');
display(GOTable(ix(1:40),:))
