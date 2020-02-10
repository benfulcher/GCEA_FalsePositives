% See if we can get more nuanced spatial scoring going
%-------------------------------------------------------------------------------
whatSpecies = 'mouse';
whatStructFilt = 'all';
params = GiveMeDefaultParams(whatSpecies,whatStructFilt);

% Get gene-expression data:
params.g.normalizationGene = 'zscore';
params.g.normalizationRegion = 'zscore';
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
numGenes = height(geneInfo);

% Get GO categories:
GOTable = GiveMeGOData(params,geneInfo.entrez_id);
numGenes = height(geneInfo);
numGOCategories = height(GOTable);

% Get pairwise distances:
distMat = GiveMeDistanceMatrix(params.humanOrMouse,params.c.structFilter);
getUpperDiag = @(x) x(triu(true(size(x)),+1));
distUpper = getUpperDiag(distMat);

%-------------------------------------------------------------------------------
% Look at spatial scale and strength of each GO category:
d0 = nan(numGOCategories,1);
A = nan(numGOCategories,1);
B = nan(numGOCategories,1);

for j = 1:numGOCategories
    fprintf(1,'Category %u/%u\n',j,numGOCategories);
    hereIam = ismember(geneInfo.entrez_id,GOTable.annotations{j});
    geneDataZ = BF_NormalizeMatrix(geneData(:,hereIam),'zscore');
    geneDataZZ = BF_NormalizeMatrix(geneDataZ','zscore')';
    G = corr(geneDataZZ','rows','pairwise');
    upperMask = triu(true(size(G)),+1);
    try
        [f_handle,Stats,c] = GiveMeFit(distUpper,getUpperDiag(G),'exp',false);
        d0(j) = 1/c.n; % spatial scale of transcriptional autocorrelation in this category
        A(j) = c.A;
        B(j) = c.B;
    end
end
