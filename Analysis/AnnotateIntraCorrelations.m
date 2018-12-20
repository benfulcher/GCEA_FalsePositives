function GOTable = AnnotateIntraCorrelations(params,GOTable,analName)
% Annotates intracorrelation scores
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs:
if nargin < 2
    GOTable = [];
end
if nargin < 3
    analName = 'intracorr';
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
% Record within-category correlation scores:
categoryScoresRaw = nan(numGOCategories,1);
categoryScoresAbs = nan(numGOCategories,1);
parfor j = 1:numGOCategories
    matchMe = ismember(geneInfo.entrez_id,GOTable.annotations{j});
    if sum(matchMe) <= 1
        continue
    end
    geneMatrix = geneData(:,matchMe);
    [categoryScoresRaw(j),categoryScoresAbs(j)] = IntraCorrelationScore(geneMatrix);
end

%-------------------------------------------------------------------------------
GOTable.(analName) = categoryScoresRaw;
GOTable.(sprintf('%s_abs',analName)) = categoryScoresAbs;

end
