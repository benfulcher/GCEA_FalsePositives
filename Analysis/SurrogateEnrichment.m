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
fakeFlags = struct();
fakeFlags.mouse = 'surrogate-mouse';
fakeFlags.human = 'surrogate-human';
params.g.humanOrMouse = fakeFlags.(whatSpecies);
params.g.whatSurrogate = whatSurrogate;
geneDataNull = LoadMeG(params.g);
if numMaps > size(geneDataNull,2)
    error('There aren''t enough null maps to compare against...');
end
if size(geneDataNull,1)~=numAreas
    error('Different parcellation???')
end
%-------------------------------------------------------------------------------
% Get a generic GO Table:
GOTableGeneric = GetFilteredGOData(params.e.dataSource,params.e.processFilter,...
                                    params.e.sizeFilter,geneInfoReal.entrez_id);
numGOCategories = height(GOTableGeneric);

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Enrichment of genes with a given random (or spatially autocorrelated) map
numGenesReal = height(geneInfoReal);
surrogatePVals = zeros(numGOCategories,numMaps);
for i = 1:numMaps
    fprintf(1,'%u/%u\n',i,numMaps);
    map_i = geneDataNull(:,i);
    geneScores = zeros(numGenesReal,1);
    for j = 1:numGenesReal
        geneScores(j) = corr(map_i,geneDataReal(:,j),'type','Spearman','rows','pairwise');
        if isnan(geneScores(j))
            keyboard
        end
    end
    % Store random-gene enrichment results:
    GOTable_i = SingleEnrichment(geneScores,geneInfoReal.entrez_id,params.e);
    % Map to the generic table:
    [~,ia,ib] = intersect(GOTableGeneric.GOID,GOTable_i.GOID,'stable');
    if ~(ia==1:height(GOTableGeneric))
        error('Could not map to generic table');
    end
    surrogatePVals(:,i) = GOTable_i.pVal(ib);
    % List GO categories with significant p-values:
    numSig = sum(GOTable_i.pValCorr < params.e.sigThresh);
    fprintf(1,'Iteration %u/%u (%s-%s): %u GO categories have p_corr < %.2f\n',...
                i,numMaps,whatSpecies,whatSurrogate,numSig,params.e.sigThresh);
    display(GOTable_i(1:numSig,:));
end

%-------------------------------------------------------------------------------
% Save out
fileNameOut = sprintf('SurrogateGOTables_%u_%s_%s.mat',numMaps,whatSpecies,whatSurrogate);
save(fileNameOut,'GOTableGeneric','surrogatePVals','-v7.3');

end
