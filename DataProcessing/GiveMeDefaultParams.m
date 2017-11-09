function params = GiveMeDefaultParams(whatParamSet,humanOrMouse)
% Idea is to give a parameter vector of defaults
%-------------------------------------------------------------------------------

if nargin < 2
    humanOrMouse = 'mouse';
    fprintf(1,'Mouse by default. Cute.\n');
end

params = struct();

switch whatParamSet
case 'conn'
    % Connectome processing
    switch humanOrMouse
    case 'mouse'
        params.connectomeSource = 'mouse-Oh'; % 'Oh-cortex'
        params.whatWeightMeasure = 'NCD';
        params.whatHemispheres = 'right';
    case 'human'
        params.connectomeSource = 'human-HCP'; % 'Oh-cortex'
        params.whatWeightMeasure = 'density';
        params.whatHemispheres = 'left';
    end
    params.pThreshold = 0.05;
    params.structFilter = 'all'; % 'cortex', 'all'

case 'gene'
    params.humanOrMouse = humanOrMouse;
    switch humanOrMouse
    case 'mouse'
        % Gene data processing
        params.energyOrDensity = 'energy'; % what gene expression data to use
        params.normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
        params.normalizationRegion = 'none'; % 'none', 'zscore'
        params.subsetOfGenes = []; %[]; only look at the first X genes.
                                % Set to empty, [], to use all genes
    case 'human'
        params.whatParcellation = 'APARC'; % 'APARC', 'HCP'
        params.probeSelection = 'variance'; % 'mean, 'variance'
        params.normalizationInternal = 'robustSigmoid'; % 'robustSigmoid', 'none'
        % Additional 'in-house' normalization:
        params.normalizationGene = 'none'; % 'none', 'mixedSigmoid'
        params.normalizationRegion = 'none'; % 'none', 'zscore'
        params.subsetOfGenes = []; %[]; only look at the first X genes.
                                % Set to empty, [], to use all genes
    end
case 'enrichment'
    % GO enrichment processing
    params.whatSource = 'direct'; % 'direct' (Direct annotations from GO)
                                  % 'GEMMA' (Annotations inferred from GEMMA)
    params.processFilter = 'biological_process';
    params.sizeFilter = [5,100];
    params.numIterations = 20000; % number of iterations for GSR
    params.enrichmentSigThresh = 0.05;

end

end
