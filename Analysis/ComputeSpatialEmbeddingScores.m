function ComputeSpatialEmbeddingScores()
% Compute spatial embedding scores (as enrichment with distance) for all GO
% categories (ready for subsequent analysis)
%-------------------------------------------------------------------------------

whatSpecies = {'mouse','mouse','human'};
whatStructFilt = {'all','cortex','cortex'};
doSave = true;

for s = 1:3
    params = GiveMeDefaultParams(whatSpecies{s},whatStructFilt{s});
    % Compute spatial autocorrelation scores per GO category:
    params.g.normalizationGene = 'zscore';
    params.g.normalizationRegion = 'zscore';
    params.e.numNullSamples = 10; % for speed since we don't actually use the p-values
    results = geneEnrichmentDistance(params,false,doSave);
end

end
