
%-------------------------------------------------------------------------------
% Connectome parameters
%-------------------------------------------------------------------------------
connectomeSource = 'Oh'; % 'Oh-cortex'
pThreshold = 0.05;
whatHemispheres = 'right';
whatWeightMeasure = 'NCD';

[A_bin,regionStruct,adjPVals] = GiveMeAdj(connectomeSource,pThreshold,true,...
                                whatWeightMeasure,whatHemispheres,justCortex);

%-------------------------------------------------------------------------------
% Gene parameters
%-------------------------------------------------------------------------------
[geneData,geneInfo,structInfo] = LoadMeG({'none','none'},'energy');
geneEntrezIDs = geneInfo.entrez_id;
gScore = nanmean(geneData,1);
gScore = gScore(randperm(length(gScore)));


%-------------------------------------------------------------------------------
% Degree enrichment:
%-------------------------------------------------------------------------------
k = sum(A_bin,1)' + sum(A_bin,2);
gScore = zeros(height(geneInfo),1);
for i = 1:height(geneInfo)
    gScore(i) = corr(k,geneData(:,i),'type','Spearman','rows','pairwise');
end
[GOTable,geneEntrezAnnotations] = SingleEnrichment(gScore,geneEntrezIDs,'biological_process',[5,200],20000);

%-------------------------------------------------------------------------------
% Cortical enrichment
%-------------------------------------------------------------------------------
gScore = zeros(height(geneInfo),1);
isCortex = strcmp(structInfo.divisionLabel,'Isocortex');
for i = 1:height(geneInfo)
    gScore(i) = ttest2(geneData(isCortex,i),geneData(~isCortex,i),'VarType','Unequal');
end
[GOTable,geneEntrezAnnotations] = SingleEnrichment(gScore,geneEntrezIDs,'biological_process',[5,200],20000);

%-------------------------------------------------------------------------------
% Get region-based anatomical nulls:
%-------------------------------------------------------------------------------
[GOTable,geneEntrezAnnotations] = GetFilteredGOData('biological_process',[5,200],geneInfo.entrez_id);
sizeGOCategories = cellfun(@length,geneEntrezAnnotations);
numGOCategories = height(GOTable);
numNulls = 100;
categoryScores = nan(numGOCategories,numNulls+1);
k = sum(A_bin,1)' + sum(A_bin,2);
parfor n = 1:numNulls+1
    fprintf(1,'Null %u/%u\n',n,numNulls+1);
    if n == 1
        permVector = 1:size(geneData,1);
    else
        permVector = randperm(size(geneData,1));
        % permVector = AnatomyShuffle(structInfo.divisionLabel);
    end

    gScore = zeros(height(geneInfo),1);
    for i = 1:height(geneInfo)
        gScore(i) = corr(k,geneData(permVector,i),'type','Pearson','rows','pairwise');
    end

    % Record mean scores for each category:
    for j = 1:numGOCategories
        matchMe = ismember(geneInfo.entrez_id,geneEntrezAnnotations{j});
        if sum(matchMe) <= 1
            continue
        end
        categoryScores(j,n) = nanmean(gScore(matchMe));
    end
end

%-------------------------------------------------------------------------------
% Estimate p-values:
%-------------------------------------------------------------------------------
[meanNull,stdNull,pValsPerm,pValsZ,pValsZ_corr] = EstimatePVals(categoryScores,numNulls,whatTail);
ListCategories(geneInfo,GOTable,geneEntrezAnnotations,meanNull,pValsZ,pValsZ_corr);
