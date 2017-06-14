function [GOTable,gScore,geneEntrezIDs] = EdgeEnrichment(whatEdgeMeasure,...
                                            onlyOnEdges,correctDistance,absType)
% Computes correlation between a pairwise measure and gene expression outer
% product
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Fixed processing parameters:
%-------------------------------------------------------------------------------

if nargin < 1
    whatEdgeMeasure = 'distance';
end
if nargin < 2
    onlyOnEdges = false;
end
if nargin < 3
    correctDistance = false; % false,true;
end
if nargin < 4
    absType = 'pos'; % 'pos','neg','abs' -> e.g., pos -> coexpression contribution increases with the statistic
end

%===============================================================================
% Connectome processing
connectomeSource = 'Oh'; % 'Oh-cortex'
pThreshold = 0.05;
whatHemispheres = 'right';
whatWeightMeasure = 'NCD';
whatFilter = 'all'; % 'cortex', 'all'

% Gene processing
energyOrDensity = 'energy'; % what gene expression data to use
normalizationGene = 'zscore'; % 'none', 'mixedSigmoid'
normalizationRegion = 'none'; % 'none', 'zscore'

% Null model:

% Correlations:
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
corrType = 'Spearman'; % {'Spearman','Pearson'};
pValOrStat = 'stat'; % 'pval','stat'

% Enrichment parameters:
numIterations = 20000; % number of iterations for GSR
enrichmentSigThresh = 0.05;

%===============================================================================
% Get data:
[A_wei,regionAcronyms,A_p] = GiveMeAdj(connectomeSource,pThreshold,false,...
                                    whatWeightMeasure,whatHemispheres,whatFilter);
A_bin = (A_wei~=0);
[geneData,geneInfo,structInfo] = LoadMeG({normalizationGene,normalizationRegion},energyOrDensity);
[A_bin,geneData,structInfo,keepInd] = filterStructures(whatFilter,structInfo,A_bin,geneData);
A_wei = A_wei(keepInd,keepInd);
A_p = A_p(keepInd,keepInd);

%-------------------------------------------------------------------------------
% Compute the edge measure:
if strcmp(whatEdgeMeasure,'distance')
    C = load('Mouse_Connectivity_Data.mat','Dist_Matrix');
    edgeData = C.Dist_Matrix{1,1};
    if onlyOnEdges
        edgeData(A_bin==0) = 0;
    else
        edgeData(tril(true(size(edgeData)))) = 0;
    end
else
    edgeData = GiveMeEdgeMeasure(whatEdgeMeasure,A_bin,A_wei,onlyOnEdges,A_p);
end

%-------------------------------------------------------------------------------
% Regress distance?:
if correctDistance
    C = load('Mouse_Connectivity_Data.mat','Dist_Matrix');
    distanceRegressor = C.Dist_Matrix{1,1};
    fprintf(1,'Regressing ipsilateral distances\n');
else
    fprintf(1,'No distance regressor used\n');
    distanceRegressor = []; % just compute normal correlations
end

%-------------------------------------------------------------------------------
% Compute gene scores:
%-------------------------------------------------------------------------------
[gScore,geneEntrezIDs] = GiveMeGCC(edgeData,geneData,...
                            geneInfo.entrez_id,corrType,distanceRegressor,absType,...
                            thresholdGoodGene,pValOrStat);

% Filter:
filterMe = isnan(gScore);
if any(filterMe)
    warning('What? NaNs? Should have been filtered within GiveMeGCC...')
    keyboard
    % fprintf(1,'%u genes had NaN scores and are being filtered out before enrichment\n',...
    %                         sum(filterMe));
    % gScore(filterMe) = [];
    % geneEntrezIDs(filterMe) = [];
end

%-------------------------------------------------------------------------------
% Enrichment using our in-house random-gene null method:
%-------------------------------------------------------------------------------
[GOTable,geneEntrezAnnotations] = SingleEnrichment(gScore,geneEntrezIDs,...
                                    'biological_process',[5,200],numIterations);

%-------------------------------------------------------------------------------
% ANALYSIS:
%-------------------------------------------------------------------------------
numSig = sum(GOTable.pVal_corr < 0.05);
fprintf(1,'%u significant categories at p_corr < 0.05\n',numSig);
display(GOTable(1:numSig,:));

end
