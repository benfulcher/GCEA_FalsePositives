function GOTable = AnnotateIntraCorrelations(params,GOTable,analName)
% Annotates intracorrelation scores
%-------------------------------------------------------------------------------

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
    GOTable = GetFilteredGOData(params.e.dataSource,params.e.processFilter,params.e.sizeFilter,...
                                    geneInfo.entrez_id);
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
    % Correlation matrix:
    C = corr(geneData(:,matchMe),'rows','pairwise','type','Pearson');
    % (consider looking also at negative correlations)
    corrVect = C(triu(true(size(C)),+1));
    categoryScoresRaw(j) = nanmean(corrVect);
    categoryScoresAbs(j) = nanmean(abs(corrVect));
end

%-------------------------------------------------------------------------------
GOTable.(analName) = categoryScores;


end
