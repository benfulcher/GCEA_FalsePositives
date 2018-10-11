% Determine enrichment signatures of different areas

%-------------------------------------------------------------------------------
% Scoring genes by differential expression between cortical and non-cortical areas
% Store results for mouse and human in a structure:
resultsTablesCortexDiff = struct();

% Mouse:
params = GiveMeDefaultParams('mouse');
resultsTablesCortexDiff.mouse = NodeSimpleEnrichment('isocortex','all','',params);

% Human-APARC (contains subcortical areas):
params = GiveMeDefaultParams('human');
params.g.whatParcellation = 'APARC';
params.g.normalizationInternal = 'none';
params.g.normalizationRegion = 'none';
params.g.normalizationGene = 'zscore';
resultsTablesCortexDiff.human_APARC = NodeSimpleEnrichment('isocortex','all','',params);

% Plot the combined table:
thresholdSig = 0.5;
PlotEnrichmentTables(resultsTablesCortexDiff,thresholdSig);
