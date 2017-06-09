% Idea is to shuffle node properties across nodes, generating a null
% distribution for each category separately
numNulls = 200;
whatEnrichment = 'degree';
whatShuffle = 'anatomyTwo'; % 'anatomy', 'full'
structureFilter = 'cortex';

%-------------------------------------------------------------------------------
% Load data:
%-------------------------------------------------------------------------------
[A_bin,regionAcronyms,adjPVals] = GiveMeAdj('Oh',0.05,true,'NCD','right');
[geneData,geneInfo,structInfo] = LoadMeG({'none','none'},'energy');
if strcmp(structureFilter,'cortex')
    keepStruct = strcmp(structInfo.divisionLabel,'Isocortex');
    geneData = geneData(keepStruct,:);
    structInfo = structInfo(keepStruct,:);
    A_bin = A_bin(keepStruct,keepStruct);
end
[GOTable,geneEntrezAnnotations] = GetFilteredGOData('biological_process',[5,100],geneInfo.entrez_id);
sizeGOCategories = cellfun(@length,geneEntrezAnnotations);
numGOCategories = height(GOTable);

%-------------------------------------------------------------------------------
% Define gene scoring function
%-------------------------------------------------------------------------------
switch whatEnrichment
case 'degree'
    % Score based on correlations to degree
    k = sum(A_bin,1)' + sum(A_bin,2);
    score_fn = @(x) corr(k,x,'type','Spearman','rows','pairwise');
end

%-------------------------------------------------------------------------------
% Define shuffle function
%-------------------------------------------------------------------------------
switch whatShuffle
case 'anatomyTwo'
    shuffle_fn = @()AnatomyShuffle(structInfo.divisionLabel,'twoBroad');
case 'anatomyFive'
    shuffle_fn = @()AnatomyShuffle(structInfo.divisionLabel,'fiveByEye');
case 'all'
    shuffle_fn = @()randperm(size(geneData,1));
end

%-------------------------------------------------------------------------------
% Assign scores to categories of genes
%-------------------------------------------------------------------------------
categoryScores = nan(numGOCategories,numNulls+1);
parfor n = 1:numNulls+1
    fprintf(1,'Null %u/%u\n',n,numNulls+1);
    if n == 1
        permVector = 1:size(geneData,1);
    else
        permVector = shuffle_fn();
    end

    gScore = zeros(height(geneInfo),1);
    for i = 1:height(geneInfo)
        gScore(i) = score_fn(geneData(permVector,i));
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
