%-------------------------------------------------------------------------------
% Fixed processing parameters:
%-------------------------------------------------------------------------------

whatEdgeMeasure = 'distance';
onlyOnEdges = true;

% Connectome processing
connectomeSource = 'Oh'; % 'Oh-cortex'
pThreshold = 0.05;
whatHemispheres = 'right';
whatWeightMeasure = 'NCD';
whatFilter = 'all'; % 'cortex', 'all'

% Gene processing
energyOrDensity = 'energy'; % what gene expression data to use
normalizationGene = 'none'; % 'none', 'mixedSigmoid'
normalizationRegion = 'none'; % 'none', 'zscore'

% Null model:

% Correlations:
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
corrType = 'Spearman'; % {'Spearman','Pearson'};
pValOrStat = 'stat'; % 'pval','stat'
absType = 'neg'; % 'pos','neg','abs' -> e.g., pos -> coexpression contribution increases with the statistic
correctDistance = false; % false,true;

% Enrichment parameters:
numIterations = 10000; % number of iterations for GSR
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
[gScore,geneEntrezIDs] = GiveMeGCC(edgeData,geneData,...
                            geneInfo.entrez_id,corrType,distanceRegressor,absType,...
                            thresholdGoodGene,pValOrStat);

% Filter first:
filterMe = isnan(gScore);
gScore(filterMe) = [];
geneEntrezIDs(filterMe) = [];

% Enrichment using our in-house random-gene null method:
[GOTable,geneEntrezAnnotations] = SingleEnrichment(gScore,geneEntrezIDs,...
                                    'biological_process',[5,200],numIterations);

% ANALYSIS:
numSig = sum(GOTable.pVal_corr < 0.05);
fprintf(1,'%u significant categories at p_corr < 0.05\n',numSig);
display(GOTable(1:numSig,:));
