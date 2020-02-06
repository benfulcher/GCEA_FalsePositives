function params = GiveMeDefaultParams(humanOrMouse,structFilter)
% Produces a parameter vector of default values for processing and analysis
%
% Gives a unified parameter structure broken into 4 components:
% * c (connectome data processing)
% * g (gene-expression data processing)
% * e (GO enrichment settings)
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

if nargin < 2
    % Default structure filter:
    switch humanOrMouse
    case {'human','surrogate-human'}
        structFilter = 'cortex'; % 'cortex', 'all'
    case {'mouse','surrogate_mouse'}
        structFilter = 'all'; % 'cortex', 'all'
    end
end
params.structFilter = structFilter;

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
    params.c.connectomeSource = 'human'; % 'human-HCP-APARC', 'human-HCP-HCP'
    params.c.whatDensity = 15; % can be 15 or 25
    params.c.whatWeightMeasure = 'density';
    params.c.whatHemispheres = 'left';
end
params.c.pThreshold = 0.05;
params.c.structFilter = params.structFilter;

%-------------------------------------------------------------------------------
% Gene expression processing options
%-------------------------------------------------------------------------------
params.g = struct();
params.g.humanOrMouse = humanOrMouse;
params.g.structFilter = params.structFilter;
params.g.whatSurrogate = 'spatialLag';
params.g.minGoodPropGene = 0.5;
params.g.minGoodPropArea = 0.5;
switch humanOrMouse
case 'mouse'
    params.g.energyOrDensity = 'energy'; % what gene expression data to use
    params.g.normalizationGene = 'none'; % 'none', 'zscore', 'mixedSigmoid'
    params.g.normalizationRegion = 'none'; % 'none', 'zscore'
    params.g.subsetOfGenes = []; % only analyze the first X genes
                                 % Set to empty, [], to use all genes
case 'human'
    % New data (Aurina 2018) does not allow selection of these options:
    % params.g.probeSelection = 'DS'; % 'mean, 'variance', 'DS'
    % params.g.normalizationInternal = 'robustSigmoid'; % 'robustSigmoid', 'none'
    params.g.whatParcellation = 'cust100'; % 'HCP', 'cust100'
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
    params.e.dataSource = 'mouse-direct';
case 'human'
    params.e.dataSource = 'human-direct';
end

params.e.processFilter = 'biological_process';
params.e.sizeFilter = [10,200];
params.e.numNullSamples = 40000; % number of null samples when computing gene-score significance
params.e.sigThresh = 0.05; % display categories with corrected p-value below this threshold
% params.e.sizeFix = []; % set the number of annotations to all categories to be a fixed value
                       % (mainly useful for )

%-------------------------------------------------------------------------------
% Parameters specific to ensemble enrichment:
params.e.whatCorr = 'Spearman';
params.e.aggregateHow = 'mean';
params.e.whatEnsemble = 'randomMap'; % 'randomMap', 'customEnsemble'
% *For custom ensemble:
params.e.useAutoSpatial = true;
% dataFileSurrogate: point to the file containing the custom maps)
% (this information is only used for 'customEnsemble')
params.e.dataFileSurrogate = FindNullFile(params);
% Filename to save results out to:
params.e.fileNameOut = GiveMeEnsembleEnrichmentOutputFileName(params);

%-------------------------------------------------------------------------------
% Properties of nulls
%-------------------------------------------------------------------------------
params.nulls.numNullsFPSR = 10000;
params.nulls.permTestP = false; % permutation test or gaussian-approximation

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
