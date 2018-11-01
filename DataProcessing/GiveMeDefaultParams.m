function params = GiveMeDefaultParams(humanOrMouse)
% Produces a parameter vector of default values for processing and analysis
%
% Gives a unified parameter structure broken into 4 components:
% * c (connectome processing)
% * g (gene expression processing)
% * e (GO enrichment processing)
% * gcc (settings for computing GCC scores)
%-------------------------------------------------------------------------------

% Check inputs:
if nargin < 1 || isempty(humanOrMouse)
    humanOrMouse = 'mouse';
    fprintf(1,'Mouse by default. Cute.\n');
end

% Initialize:
params = struct();
params.humanOrMouse = humanOrMouse;

% Filter structures:
structFilter = 'all'; % 'isocortex', 'all'

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
params.c.structFilter = structFilter;

%-------------------------------------------------------------------------------
% Gene expression processing options
%-------------------------------------------------------------------------------
params.g = struct();
params.g.humanOrMouse = humanOrMouse;
params.g.useSurrogate = false;
params.g.structFilter = structFilter;
switch humanOrMouse
case 'mouse'
    params.g.energyOrDensity = 'energy'; % what gene expression data to use
    params.g.normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
    params.g.normalizationRegion = 'zscore'; % 'none', 'zscore'
    params.g.subsetOfGenes = []; % only analyze the first X genes
                                 % Set to empty, [], to use all genes
case 'human'
    % New data (Aurina 2018) does not allow selection of these options:
    % params.g.whatParcellation = 'cust100'; % 'APARC', 'HCP', 'cust100'
    % params.g.probeSelection = 'DS'; % 'mean, 'variance', 'DS'
    % params.g.normalizationInternal = 'robustSigmoid'; % 'robustSigmoid', 'none'
    params.g.normalizeSeparately = true; % whether to normalize cortex/subcortex separately

    % Additional 'in-house' normalization:
    params.g.normalizationGene = 'none'; % 'none', 'mixedSigmoid'
    params.g.normalizationRegion = 'none'; % 'none', 'zscore'
    params.g.subsetOfGenes = []; % only analyze the first X genes
                                 % Set to empty, [], to use all genes
end

%-------------------------------------------------------------------------------
% GO enrichment options
%-------------------------------------------------------------------------------
params.e = struct();
params.e.humanOrMouse = humanOrMouse;

% 1) Specify a data source:
switch humanOrMouse
case 'mouse'
    params.e.dataSource = 'mouse-direct'; % 'mouse-direct', 'mouse-GEMMA'
case 'human'
    params.e.dataSource = 'human-direct';
end

params.e.processFilter = 'biological_process';
params.e.sizeFilter = [10,100];
params.e.numSamples = 20000; % number of null samples when computing gene score significance
params.e.sigThresh = 0.05; % display categories with corrected p-value below this threshold

%-------------------------------------------------------------------------------
% Computing GCC scores
%-------------------------------------------------------------------------------
params.gcc = struct();
params.gcc.onlyConnections = false; % only look where there are structural connections
params.gcc.regressDistance = false; % whether to regress distance
% Settings for computing correlations:
params.gcc.whatCorr = 'Spearman'; % 'Pearson', 'Spearman'
params.gcc.pValOrStat = 'stat'; % 'pval','stat'
params.gcc.thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
params.gcc.absType = 'neg'; % 'pos','neg','abs' -> e.g., pos -> coexpression contribution increases with the statistic

% Computing nulls through shuffling:
params.gcc.numNulls = 50; % number of nulls
params.gcc.whatTail = 'right'; % right-tailed p-values


end
