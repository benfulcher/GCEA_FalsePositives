

% Run enrichment:
params = GiveMeDefaultParams('mouse');
% To make GCC scores make sense -- expression needs to be [0,1] normalized:
params.g.normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
params.g.normalizationRegion = 'zscore'; % 'none', 'zscore'
params.c.structFilter = 'all';

%-------------------------------------------------------------------------------
% Enrichment of random maps with distance
% (nothing, as expected)
params.g.humanOrMouse = 'surrogate-mouse';
mouse_surrogate_enrichment = geneEnrichmentDistance(params);

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Enrichment of genes with a given random (spatially correlated) map
numMaps = 3;
params.g.humanOrMouse = 'mouse';
[geneDataReal,geneInfoReal,structInfoReal] = LoadMeG(params.g);
params.g.humanOrMouse = 'surrogate-mouse';
[geneDataNull,geneInfoNull,structInfoNull] = LoadMeG(params.g);
% Score each gene with map 1:
numGenesReal = height(geneInfoReal);
GOTables = cell(numMaps,1);
for i = 1:numMaps
    map_i = geneDataNull(:,i);
    geneScores = zeros(numGenesReal,1);
    for j = 1:numGenesReal
        geneScores(j) = corr(map_i,geneDataReal(:,j),'type','Spearman');
    end
    GOTables{i} = SingleEnrichment(geneScores,geneInfoReal.entrez_id,params.e);
    %-------------------------------------------------------------------------------
    % List GO categories with significant p-values:
    numSig = sum(GOTable.pValCorr < params.e.sigThresh);
    fprintf(1,'%u GO categories have p_corr < %.2f\n',numSig,params.e.sigThresh);
    display(GOTable(1:numSig,:));
end

%-------------------------------------------------------------------------------
% Save out
fileNameOut = 'SurrogateGOTables.mat';
save(fileNameOut,GOTables);
