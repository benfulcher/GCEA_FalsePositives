function params = GiveMeDefaultParams(whatParamSet)
% Idea is to give a parameter vector of defaults
%-------------------------------------------------------------------------------

params = struct();

switch whatParamSet
case 'conn'
    % Connectome processing
    params.connectomeSource = 'Oh'; % 'Oh-cortex'
    params.pThreshold = 0.05;
    params.whatHemispheres = 'right';
    params.whatWeightMeasure = 'NCD';
    params.structFilter = 'all'; % 'cortex', 'all'

case 'gene'
    % Gene data processing
    params.energyOrDensity = 'energy'; % what gene expression data to use
    params.normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
    params.normalizationRegion = 'zscore'; % 'none', 'zscore'
    params.subsetOfGenes = []; %[]; only look at the first X genes.
                                % Set to empty, [], to use all genes

case 'enrichment'
    % GO enrichment processing
    params.processFilter = 'biological_process';
    params.sizeFilter = [5,100];
    params.numIterations = 20000; % number of iterations for GSR
    params.enrichmentSigThresh = 0.05;

end

end
