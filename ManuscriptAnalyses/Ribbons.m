% Ribbons
% Put multiple ribbons together
%-------------------------------------------------------------------------------

% Store results tables in this struct:
results = struct();

%===============================================================================
% Distance
%===============================================================================
% Mouse brain:
params = GiveMeDefaultParams('mouse');
params.e.numSamples = 100000;
params.g.normalizationGene = 'zscore';
params.g.normalizationRegion = 'zscore';
results.mouseDistance = geneEnrichmentDistance(params);

% Human cortex:
params = GiveMeDefaultParams('human');
params.e.numSamples = 100000;
params.g.normalizationGene = 'zscore';
params.g.normalizationRegion = 'zscore';
results.humanDistance = geneEnrichmentDistance(params);

%===============================================================================
% Within-category correlation
%===============================================================================
numSamples = 500;
params = GiveMeDefaultParams('mouse');
results.mouseIntra = IntraCorrelationByCategory(params,'geneShuffle',numSamples);
results.mouseIntra.pValCorr = results.mouseIntra.pValZCorr;
params = GiveMeDefaultParams('human');
results.humanIntra = IntraCorrelationByCategory(params,'geneShuffle',numSamples);
results.humanIntra.pValCorr = results.humanIntra.pValZCorr;

%-------------------------------------------------------------------------------
% Are coexpression scores related to spatial correlation scores?
justMouseResults = struct('mouseDistance',results.mouseDistance,'mouseIntra',results.mouseIntra);
% e.g., for mouse:
[rowVectorResults,GOTerms,allGOIDs,tableNames] = CombineTables(justMouseResults,'mouse',...
    {'pValZ','pValZ'});
% {'meanScore','mouse'});
plot(rowVectorResults(1,:),rowVectorResults(2,:),'.k')
xlabel(tableNames{1})
ylabel(tableNames{2})

%===============================================================================
% Significant categories under the spatial lag model
%===============================================================================
results.mouseSpatialLag = SurrogateEnrichmentProcess('mouse',1000,'spatialLag');
results.humanSpatialLag = SurrogateEnrichmentProcess('human','spatialLag');

%===============================================================================
% Can we explain the final result in terms of the former characteristics?
[rowVectorResults,GOTerms,allGOIDs,allTableNames] = CombineTables(results,'mouse',...
    {'pValZ','meanScore','pValZ','human_abs','sumUnderSig','sumUnderSigx'});
f = figure('color','w');
isGood = (sum(~isfinite(rowVectorResults),1)==0);
scatter(rowVectorResults(1,isGood),...
        rowVectorResults(3,isGood),...
        rowVectorResults(5,isGood)+1,...
        rowVectorResults(5,isGood),'filled');


mouseScoreDistance = -log10(rowVectorResults(1,isGood));
mouseScoreDistance(isinf(mouseScoreDistance)) = 6;
mouseScoreCorr = -log10(rowVectorResults(3,isGood));
mouseScoreCorr(isinf(mouseScoreCorr)) = 6;

scatter(mouseScoreDistance,...
        mouseScoreCorr,...
        rowVectorResults(5,isGood)+1,...
        rowVectorResults(5,isGood),'filled');
colormap('jet')

%-------------------------------------------------------------------------------
% Check em out together:
thresholdSig = [0.05,2,2];
[allGOIDsSort,GONamesSort] = PlotEnrichmentTables(results,thresholdSig);
