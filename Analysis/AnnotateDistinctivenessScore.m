function GOTable = AnnotateDistinctivenessScore(params,GOTable,analName)
% Annotates distinctiveness scores
%-------------------------------------------------------------------------------
% Categories have high distinctiveness that have transcriptional maps that are
% less related to those of other genes in the brain
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs:
if nargin < 2
    GOTable = [];
end
if nargin < 3
    analName = 'distinctiveness';
end

%-------------------------------------------------------------------------------
% Load expression data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);

%-------------------------------------------------------------------------------
% Load annotation data:
if isempty(GOTable)
    GOTable = GiveMeGOData(params,geneInfo.entrez_id);
end
numGOCategories = height(GOTable);

%-------------------------------------------------------------------------------
% Compute coexpression between all genes:
coexp = corr(geneData,'type','Pearson','rows','pairwise');

%-------------------------------------------------------------------------------
% Record within-category correlation scores:
categoryScoresDistinctRaw = nan(numGOCategories,1);
categoryScoresDistinctAbs = nan(numGOCategories,1);
for j = 1:numGOCategories
    matchMe = ismember(geneInfo.entrez_id,GOTable.annotations{j});
    if sum(matchMe) <= 1
        continue
    end
    % Correlation to genes outside the category
    coexpExt = coexp(matchMe,~matchMe);
    categoryScoresDistinctRaw(j) = nanmean(coexpExt(:));
    categoryScoresDistinctAbs(j) = nanmean(abs(coexpExt(:)));
end

%-------------------------------------------------------------------------------
GOTable.(sprintf('%s_raw',analName)) = categoryScoresDistinctRaw;
GOTable.(sprintf('%s_abs',analName)) = categoryScoresDistinctAbs;

end
