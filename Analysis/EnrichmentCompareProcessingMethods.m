% Idea is to compare different preprocessing steps on the enrichment of a given
% edge property
whatEdgeProperty = 'betweenness';


%-------------------------------------------------------------------------------
% Fixed parameters:
energyOrDensity = 'energy'; % what gene expression data to use
pValOrStat = 'stat'; % 'pval','stat'
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
numIterationsErmineJ = 20000; % number of iterations for GSR in ermineJ

%-------------------------------------------------------------------------------
% Set up a structure array containing all of the different processing options:
absTypes = {true,false};
corrTypes = {'Spearman','Pearson'};
normalizationGeneTypes = {'none','log10','robustSigmoid'};
normalizationRegionTypes = {'none','zscore','robustSigmoid'};
correctDistanceTypes = {true,false};

cntr = 0;
for i = 1:length(absTypes)
    absType = absTypes(i);
    for j = 1:length(corrTypes)
        corrType = corrTypes{j};
        for k = 1:length(normalizationGeneTypes)
            normalizationGeneType = normalizationGeneTypes{k};
            for l = 1:length(normalizationRegionTypes)
                normalizationRegionType = normalizationRegionTypes{l};
                for m = 1:length(correctDistanceTypes)
                    correctDistance = correctDistanceTypes{m};
                    %----------------------------------------------------------
                    cntr = cntr + 1;
                    processingSteps(cntr).abs = absType;
                    processingSteps(cntr).corrType = corrType;
                    processingSteps(cntr).normalizationGene = normalizationGeneType;
                    processingSteps(cntr).normalizationRegion = normalizationRegionType;
                    processingSteps(cntr).correctDistance = correctDistance;
                end
            end
        end
    end
end
numProcessingTypes = length(processingSteps);
fprintf(1,'Comparing %u different processing parameters\n');

%-------------------------------------------------------------------------------
% Get edge data:
%-------------------------------------------------------------------------------
C = load('Mouse_Connectivity_Data.mat'); % C stands for connectome data
pThreshold = 0.05;
A_bin = GiveMeAdj(C,'binary','ipsi',0,pThreshold);
edgeData = edge_betweenness_bin(A_bin);
f = figure('color','w');
histogram(edgeData(edgeData~=0));
xlabel('Edge betweenness (binary)')

%-------------------------------------------------------------------------------
% Get scores for genes:
%-------------------------------------------------------------------------------
enrichmentTables = cell(numProcessingTypes);
for i = 1:numProcessingTypes
    % Load in our gene data:
    [GeneStruct,geneData] = LoadMeG(true,{processingSteps(i).normalizationGene,...
                            processingSteps(i).normalizationRegion},energyOrDensity);
    geneEntrezIDs = [GeneStruct.gene_entrez_id];
    % Compute gene scores:
    % (sometimes entrez IDs change -- e.g., when matching to distance results)
    [gScore,geneEntrezIDs] = GiveMeGCC(edgeData,geneData,geneEntrezIDs,processingSteps(i).corrType,...
                    processingSteps(i).correctDistance,thresholdGoodGene,pValOrStat);
    % Do enrichment:
    fileNameWrite = writeErmineJFile('tmp',gScore,geneEntrezIDs,'bin_betw');
    ermineJResults = RunErmineJ(fileNameWrite,numIterationsErmineJ);
    % <<<Keep ermineJResults for p-values under 0.1>>>
    enrichmentTables{i} = ermineJResults(ermineJResults.pVal_corr < 0.1);
end
