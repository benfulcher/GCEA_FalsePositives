function ComputeSpatialEmbeddingScores()
% Compute spatial embedding scores (as enrichment with distance) for all GO
% categories (ready for subsequent analysis)
%-------------------------------------------------------------------------------

whatSpecies = {'mouse','mouse','human'};
whatStructFilt = {'all','cortex','cortex'};
doSave = true;

% Compute spatial autocorrelation scores per GO category:
for s = 1:3
    params = GiveMeDefaultParams(whatSpecies{s},whatStructFilt{s});

    % Get gene-expression data:
    params.g.normalizationGene = 'zscore';
    params.g.normalizationRegion = 'zscore';
    [geneData,geneInfo,structInfo] = LoadMeG(params.g);
    numGenes = height(geneInfo);

    % Get pairwise distances:
    distMat = GiveMeDistanceMatrix(params.humanOrMouse,params.c.structFilter);
    getUpperDiag = @(x) x(triu(true(size(x)),+1));
    distUpper = getUpperDiag(distMat);

    % Compute spatial embedding scores gene-by-gene
    geneScores = zeros(numGenes,1);
    parfor i = 1:numGenes
        g = geneData(:,i);
        GCC = g*g'; % self product at each edge
        % Does the self-correlation of expression depend on distance?:
        geneScores(i) = -corr(distUpper,getUpperDiag(GCC),...
                                    'type',params.gcc.whatCorr,...
                                    'rows','pairwise');
    end

    % Compute the mean score in each GO category:
    fprintf(1,'Mean the score in each GO category!\n');
    params.e.numNullSamples = 0;
    GOTable = SingleEnrichment(geneScores,geneInfo.entrez_id,params.e);

    % -----------------------
    % Save out to .mat file:
    fileNameMat = fullfile('DataOutputs',GiveMeDistanceScoreFileName(params));
    save(fileNameMat,'GOTable','params'); % ,'geneEntrez','geneDistanceScores'
    fprintf(1,'Saved distance enrichment scores to %s\n',fileNameMat);
end

end
