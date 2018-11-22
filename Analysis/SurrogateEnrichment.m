function SurrogateEnrichment(whatSpecies,numMaps,whatSurrogate)
% Get enrichment results across surrogate spatial maps
%-------------------------------------------------------------------------------

if nargin < 1
    whatSpecies = 'mouse';
end
if nargin < 2
    numMaps = 1000; % number of null maps to test
end
if nargin < 3
    whatSurrogate = 'spatialLag';
end

%-------------------------------------------------------------------------------
% Get real data:
params = GiveMeDefaultParams(whatSpecies);
params.g.humanOrMouse = whatSpecies;
[geneDataReal,geneInfoReal,structInfoReal] = LoadMeG(params.g);
numGenes = height(geneInfoReal);
numAreas = height(structInfoReal);

%-------------------------------------------------------------------------------
% Get surrogate data, geneDataNull (each column is a null spatial map)
switch whatSurrogate
case 'spatialLag'
    % Surrogate maps pre-generated using the spatial lag model:
    fakeFlags = struct();
    fakeFlags.mouse = 'surrogate-mouse';
    fakeFlags.human = 'surrogate-human';
    params.g.humanOrMouse = fakeFlags.(whatSpecies);
    geneDataNull = LoadMeG(params.g);
case 'randomShuffle'
    % Surrogate maps generated through random shuffling across brain areas:
    geneDataNull = geneDataReal;
    for j = 1:numGenes
        rp = randperm(numAreas);
        geneDataNull(:,j) = geneDataReal(rp,j);
    end
otherwise
    error('Unknown surrogate method: ''%s''',whatSurrogate);
end
if numMaps > size(geneDataNull,2)
    error('There aren''t enough null maps to compare against...');
end

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Enrichment of genes with a given random (or spatially autocorrelated) map
numGenesReal = height(geneInfoReal);
GOTables = cell(numMaps,1);
for i = 1:numMaps
    fprintf(1,'%u/%u\n',i,numMaps);
    map_i = geneDataNull(:,i);
    geneScores = zeros(numGenesReal,1);
    for j = 1:numGenesReal
        geneScores(j) = corr(map_i,geneDataReal(:,j),'type','Spearman');
    end
    % Store random-gene enrichment results:
    GOTables{i} = SingleEnrichment(geneScores,geneInfoReal.entrez_id,params.e);
    % List GO categories with significant p-values:
    numSig = sum(GOTables{i}.pValCorr < params.e.sigThresh);
    fprintf(1,'Iteration %u/%u (%s-%s): %u GO categories have p_corr < %.2f\n',
                i,numMaps,whatSpecies,whatSurrogate,numSig,params.e.sigThresh);
    display(GOTables{i}(1:numSig,:));
end

%-------------------------------------------------------------------------------
% Save out
fileNameOut = sprintf('SurrogateGOTables_%u_%s_%s.mat',numMaps,whatSpecies,whatSurrogate);
save(fileNameOut,'GOTables','-v7.3');

end
