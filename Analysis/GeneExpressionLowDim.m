function pcScore = GeneExpressionLowDim(structureFilter,normHow)
% Plot brain regions in a low-dimensional gene expression space
%-------------------------------------------------------------------------------
if nargin < 1
    structureFilter = 'cortex'; % 'cortex','all'
end
if nargin < 2
    normHow = 'zscore'; % 'zscore','scaledSigmoid'
end

%-------------------------------------------------------------------------------
% Load expression energy:
%-------------------------------------------------------------------------------
[geneData,geneInfo,structInfo] = LoadMeG({'none','none'},'energy');
if strcmp(structureFilter,'cortex')
    keepStruct = strcmp(structInfo.divisionLabel,'Isocortex');
    geneData = geneData(keepStruct,:);
    structInfo = structInfo(keepStruct,:);
end
% Normalize by scaled sigmoid:
geneDataNorm = BF_NormalizeMatrix(geneData,normHow);

%-------------------------------------------------------------------------------
% Compute Principal Components using als algorithm to deal with missing data:
%-------------------------------------------------------------------------------
fprintf(1,'Computing PCs for %ux%u matrix...',size(geneDataNorm,1),size(geneDataNorm,2));
[pcCoeff,pcScore,~,~,~] = pca(geneDataNorm,'NumComponents',2,'algorithm','als');
fprintf(1,'\n');

%-------------------------------------------------------------------------------
% Plot as a scatter:
%-------------------------------------------------------------------------------
RegionScatterPlot(structInfo,pcScore(:,1),pcScore(:,2),'geneExp-PC1','geneExp-PC2','Pearson',true);

%-------------------------------------------------------------------------------
% Enrichment? -> NodeEnrichment.m
%-------------------------------------------------------------------------------
% Can do enrichment of pc scores, implemented in NodeEnrichment


end
