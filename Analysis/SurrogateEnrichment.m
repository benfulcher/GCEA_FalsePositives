function SurrogateEnrichment(params)
% Compute conventional enrichment results across surrogate spatial maps
% (for FPSR results)
%-------------------------------------------------------------------------------
if nargin < 1
    params = GiveMeDefaultParams('mouse');
end
% params.g.whatSurrogate; % the ensemble of maps to use
numMaps = params.nulls.numNullsFPSR; % number of null maps to test against

%-------------------------------------------------------------------------------
% Get real data:
[geneDataReal,geneInfoReal,structInfoReal] = LoadMeG(params.g);
numGenes = height(geneInfoReal);
numAreas = height(structInfoReal);

%-------------------------------------------------------------------------------
% Get surrogate data, geneDataNull (each column is a null spatial map)
params.g.humanOrMouse = sprintf('surrogate-%s',params.humanOrMouse);
geneDataNull = LoadMeG(params.g);
if numMaps > size(geneDataNull,2)
    error('There aren''t enough null maps to compare against...');
end
if size(geneDataNull,1)~=numAreas
    error('Size not matching---different parcellation???')
end

%-------------------------------------------------------------------------------
% Get a generic GO Table:
GOTableGeneric = GiveMeGOData(params,geneInfoReal.entrez_id);
numGOCategories = height(GOTableGeneric);

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Enrichment of genes with a given null spatial map
numGenesReal = height(geneInfoReal);
surrogatePValsPerm = zeros(numGOCategories,numMaps);
surrogatePValsZ = zeros(numGOCategories,numMaps);
for i = 1:numMaps
    fprintf(1,'%u/%u\n',i,numMaps);
    map_i = geneDataNull(:,i);

    % Get randomized 'real' data for null distribution:
    if ~isempty(params.nulls.customShuffle)
        switch params.nulls.customShuffle
        case 'coordinatedSpatialShuffle'
            fprintf(1,'Coordinated spatial shuffle of gene-expression data at %u/%u\n',i,numMaps);
            geneDataReal = ShuffleMyMatrix(geneDataReal,'coordinatedRowShuffle');
        case 'independentSpatialShuffle'
            fprintf(1,'Independent spatial shuffle of gene-expression data at %u/%u\n',i,numMaps);
            geneDataReal = ShuffleMyMatrix(geneDataReal,'randomUniform');
        case 'none'
            % Do nothing
        otherwise
            error('Unknown custom shuffling option ''%s''',params.nulls.customShuffle);
        end
    end

    % Compute gene scores as (spatial) correlations:
    geneScores = zeros(numGenesReal,1);
    for j = 1:numGenesReal
        geneScores(j) = corr(map_i,geneDataReal(:,j),'type',params.e.whatCorr,'rows','pairwise');
        if isnan(geneScores(j))
            error('Error computing gene score for map %u and gene %u',i,j);
        end
    end

    % Store GSEA results (from conventional random-gene nulls):
    GOTable_i = SingleEnrichment(geneScores,geneInfoReal.entrez_id,params.e);

    % Map results onto a generic GO table:
    [~,ia,ib] = intersect(GOTableGeneric.GOID,GOTable_i.GOID,'stable');
    if ~(ia==1:height(GOTableGeneric))
        error('Could not map to generic table');
    end

    % Store RAW, UNCORRECTED p-values:
    % As a permutation test:
    surrogatePValsPerm(:,i) = GOTable_i.pValPerm(ib);
    % Inferred from a Gaussian approximation to the null:
    surrogatePValsZ(:,i) = GOTable_i.pValZ(ib);

    % List GO categories with significant (corrected) p-values:
    numSig = sum(GOTable_i.pValZCorr < params.e.sigThresh);
    fprintf(1,'Iteration %u/%u (%s-%s): %u GO categories have pZ_corr < %.2f\n',...
                i,numMaps,whatSpecies,whatSurrogate,numSig,params.e.sigThresh);
    display(GOTable_i(1:numSig,:));
end

%-------------------------------------------------------------------------------
% Save out to a .mat file
fileNameOut = GiveMeFPSRFileName(params)
save(fileNameOut,'GOTableGeneric','surrogatePValsPerm','surrogatePValsZ','-v7.3');
fprintf(1,'Results of %u iterations saved to %s\n',numMaps,fileNameOut);

end
