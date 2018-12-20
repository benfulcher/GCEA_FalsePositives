%-------------------------------------------------------------------------------
% Ribbons
%-------------------------------------------------------------------------------
% Put multiple ribbons together
%-------------------------------------------------------------------------------

% Store results tables in this struct:
results = struct();

%===============================================================================
% Distance
%===============================================================================
% Mouse brain:
params = GiveMeDefaultParams('mouse');
params.g.normalizationGene = 'zscore';
params.g.normalizationRegion = 'zscore';
results.mouseDistance = geneEnrichmentDistance(params);

% Human cortex:
params = GiveMeDefaultParams('human');
params.g.normalizationGene = 'zscore';
params.g.normalizationRegion = 'zscore';
results.humanDistance = geneEnrichmentDistance(params);

%===============================================================================
% Within-category correlation
%===============================================================================
numSamples = 100;
params = GiveMeDefaultParams('mouse');
results.mouseIntra = IntraCorrelationByCategory(params,'geneShuffle',numSamples);
results.mouseIntra.pValCorr = results.mouseIntra.pValZCorr;
params = GiveMeDefaultParams('human');
results.humanIntra = IntraCorrelationByCategory(params,'geneShuffle',numSamples);
results.humanIntra.pValCorr = results.humanIntra.pValZCorr;

%===============================================================================
% Significant categories under the spatial lag model
%===============================================================================
results.mouseSpatialLag = SurrogateEnrichmentProcess('mouse','spatialLag');
results.humanSpatialLag = SurrogateEnrichmentProcess('human','spatialLag');

%-------------------------------------------------------------------------------
% Check em out together:
thresholdSig = [0.05,2,2];
[allGOIDsSort,GONamesSort] = PlotEnrichmentTables(results,thresholdSig);
