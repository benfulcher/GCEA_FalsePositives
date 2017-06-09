% Idea is to shuffle node properties across nodes, generating a null
% distribution for each category separately
numNulls = 200;
whatEnrichment = 'degree';

%-------------------------------------------------------------------------------
% Get region-based anatomical nulls:
%-------------------------------------------------------------------------------
[GOTable,geneEntrezAnnotations] = GetFilteredGOData('biological_process',[5,100],geneInfo.entrez_id);
sizeGOCategories = cellfun(@length,geneEntrezAnnotations);
numGOCategories = height(GOTable);

categoryScores = nan(numGOCategories,numNulls+1);
k = sum(A_bin,1)' + sum(A_bin,2);
parfor n = 1:numNulls+1
    fprintf(1,'Null %u/%u\n',n,numNulls+1);
    if n == 1
        permVector = 1:size(geneData,1);
    else
        % permVector = randperm(size(geneData,1));
        permVector = AnatomyShuffle(structInfo.divisionLabel,'twoBroad');
    end

    gScore = zeros(height(geneInfo),1);
    for i = 1:height(geneInfo)
        gScore(i) = corr(k,geneData(permVector,i),'type','Spearman','rows','pairwise');
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
whatTail = 'right';
[meanNull,stdNull,pValsPerm,pValsZ,pValsZ_corr] = EstimatePVals(categoryScores,numNulls,whatTail);
ix_GO = ListCategories(geneInfo,GOTable,geneEntrezAnnotations,meanNull,pValsZ,pValsZ_corr);
NullSummaryPlots(pValsZ,pValsZ_corr,categoryScores,meanNull,stdNull,sizeGOCategories);
SpecificNullPlots(categoryScores,GOTable,sizeGOCategories,pValsZ_corr,ix_GO);
