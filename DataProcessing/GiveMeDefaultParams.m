function params = GiveMeDefaultParams(humanOrMouse)
%-------------------------------------------------------------------------------
% Idea is to give a parameter vector of defaults
% Gives a unified parameter structure broken into 3 components:
% * c (connectome processing)
% * g (gene expression processing)
% * e (GO enrichment processing)
%-------------------------------------------------------------------------------

% Check inputs:
if nargin < 1 || isempty(humanOrMouse)
    humanOrMouse = 'mouse';
    fprintf(1,'Mouse by default. Cute.\n');
end

% Initialize:
params = struct();
params.humanOrMouse = humanOrMouse;

%-------------------------------------------------------------------------------
% Connectome processing options
%-------------------------------------------------------------------------------
params.c = struct();
params.c.humanOrMouse = humanOrMouse;
switch humanOrMouse
case 'mouse'
    params.c.connectomeSource = 'mouse-Oh'; % 'Oh-cortex'
    params.c.whatWeightMeasure = 'NCD';
    params.c.whatHemispheres = 'right';
case 'human'
    params.c.connectomeSource = 'human-HCP-HCP'; % 'human-HCP-APARC', 'human-HCP-HCP'
    params.c.whatWeightMeasure = 'density';
    params.c.whatHemispheres = 'left';
end
params.c.pThreshold = 0.05;
params.c.structFilter = 'all'; % 'isocortex', 'all'

%-------------------------------------------------------------------------------
% Gene expression processing options
%-------------------------------------------------------------------------------
params.g = struct();
params.g.humanOrMouse = humanOrMouse;
switch humanOrMouse
case 'mouse'
    % Gene data processing
    params.g.energyOrDensity = 'energy'; % what gene expression data to use
    params.g.normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
    params.g.normalizationRegion = 'none'; % 'none', 'zscore'
    params.g.subsetOfGenes = []; %[]; only look at the first X genes.
                            % Set to empty, [], to use all genes
case 'human'
    params.g.whatParcellation = 'HCP'; % 'APARC', 'HCP'
    params.g.probeSelection = 'variance'; % 'mean, 'variance'
    params.g.normalizationInternal = 'robustSigmoid'; % 'robustSigmoid', 'none'
    % Additional 'in-house' normalization:
    params.g.normalizationGene = 'none'; % 'none', 'mixedSigmoid'
    params.g.normalizationRegion = 'none'; % 'none', 'zscore'
    params.g.subsetOfGenes = []; %[]; only look at the first X genes.
                            % Set to empty, [], to use all genes
end

%-------------------------------------------------------------------------------
% GO enrichment options
%-------------------------------------------------------------------------------
params.e = struct();
params.e.humanOrMouse = humanOrMouse;

% 1) Specify a data source: % 'direct' (Direct annotations from GO)
                            % 'GEMMA' (Annotations inferred from GEMMA)
switch humanOrMouse
case 'mouse'
    params.e.whatSource = 'mouse-direct'; % 'mouse-direct', 'mouse-GEMMA'
case 'human'
    params.e.whatSource = 'human-direct';
end

params.e.processFilter = 'biological_process';
params.e.sizeFilter = [5,200];
params.e.numIterations = 20000; % number of iterations for GSR
params.e.enrichmentSigThresh = 0.05;

end
