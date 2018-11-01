

% Main parameter:
whatSpecies = 'mouse';
numMaps = 10; % number of null maps to test

% Run enrichment:
params = GiveMeDefaultParams(whatSpecies);
switch whatSpecies
case 'mouse'
    realAndFake = {'mouse','surrogate-mouse'};
case 'human'
    realAndFake = {'human','surrogate-human'};
end

% % To make GCC scores make sense -- expression needs to be [0,1] normalized:
% params.g.normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
% params.g.normalizationRegion = 'zscore'; % 'none', 'zscore'
% params.c.structFilter = 'all';
%
% %-------------------------------------------------------------------------------
% % Enrichment of random maps with distance
% % (nothing, as expected)
% params.g.humanOrMouse = 'surrogate-mouse';
% mouse_surrogate_enrichment = geneEnrichmentDistance(params);

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Enrichment of genes with a given random (spatially correlated) map
params.g.humanOrMouse = realAndFake{1};
[geneDataReal,geneInfoReal,structInfoReal] = LoadMeG(params.g);
params.g.humanOrMouse = realAndFake{2};
[geneDataNull,geneInfoNull,structInfoNull] = LoadMeG(params.g);
if numMaps > height(geneInfoNull)
    error('There aren''t enough null maps to compare against...');
end
numGenesReal = height(geneInfoReal);
GOTables = cell(numMaps,1);
for i = 1:numMaps
    fprintf(1,'%u/%u\n',i,numMaps);
    map_i = geneDataNull(:,i);
    geneScores = zeros(numGenesReal,1);
    for j = 1:numGenesReal
        geneScores(j) = corr(map_i,geneDataReal(:,j),'type','Spearman');
    end
    GOTables{i} = SingleEnrichment(geneScores,geneInfoReal.entrez_id,params.e);
    %-------------------------------------------------------------------------------
    % List GO categories with significant p-values:
    numSig = sum(GOTables{i}.pValCorr < params.e.sigThresh);
    fprintf(1,'%u GO categories have p_corr < %.2f\n',numSig,params.e.sigThresh);
    display(GOTables{i}(1:numSig,:));
end

%-------------------------------------------------------------------------------
% Save out
fileNameOut = 'SurrogateGOTables.mat';
save(fileNameOut,'GOTables');
