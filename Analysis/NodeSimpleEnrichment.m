function [GOTable,gScore,geneEntrezAnnotations] = NodeSimpleEnrichment(enrichWhat,structFilter)
%
% ---INPUTS:
% enrichWhat = 'meanExpression'; % raw mean expression level
% enrichWhat = 'varExpression'; % raw variance of expression levels
% enrichWhat = 'cortex'; % difference in expression between cerebral cortex/others
% enrichWhat = 'genePC'; % correlation with a PC of gene expression
% enrichWhat = 'degree'; % correlation with structural connectivity degree across regions
%
% structFilter = 'cortex';
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Name an analysis, and this script loads the data and does the enrichment
%-------------------------------------------------------------------------------
if nargin < 1 || isempty(enrichWhat)
    fprintf(1,'Mean expression by default\n');
    enrichWhat = 'meanExpression'; % raw mean expression level
end
if nargin < 2 || isempty(structFilter)
    fprintf(1,'No filter: all brain regions included\n');
    structFilter = 'all';
end
doRandomize = false;

%-------------------------------------------------------------------------------
% Connectome parameters
%-------------------------------------------------------------------------------
% Get default parameter sets:
cParam = GiveMeDefaultParams('conn');
eParam = GiveMeDefaultParams('enrichment');
gParam = GiveMeDefaultParams('gene');
% Ensure no normalization:
gParam.normalizationGene = 'none';
gParam.normalizationRegion = 'none';


% Binary connectome data:
[A_bin,regionAcronyms,adjPVals] = GiveMeAdj(cParam.connectomeSource,cParam.pThreshold,true,...
                                cParam.whatWeightMeasure,cParam.whatHemispheres,cParam.structFilter);
% Gene data:
[geneData,geneInfo,structInfo] = LoadMeG(gParam);

% Filter structures:
[A_bin,geneData,structInfo] = filterStructures(structFilter,structInfo,A_bin,geneData);
numStructs = height(structInfo);
numGenes = height(geneInfo);

%===============================================================================
%% -------------------------------SCORE THE GENES-------------------------------
%===============================================================================
switch enrichWhat
case 'meanExpression'
    %-------------------------------------------------------------------------------
    % Mean expression enrichment: genes scored for mean expression across the brain
    %-------------------------------------------------------------------------------
    gScore = nanmean(geneData,1);

case 'varExpression'
    %-------------------------------------------------------------------------------
    % Expression variance enrichment: genes scored for variance of expression across brain
    %-------------------------------------------------------------------------------
    gScore = nanstd(geneData,1);

case 'degree'
    %-------------------------------------------------------------------------------
    % Degree enrichment: genes scored for number of connections they make to other regions
    %-------------------------------------------------------------------------------
    k = sum(A_bin,1)' + sum(A_bin,2);
    gScore = zeros(numGenes,1);
    for i = 1:numGenes
        gScore(i) = corr(k,geneData(:,i),'type','Spearman','rows','pairwise');
    end
    % Plot:
    [~,iy] = sort(gScore,'descend');
    [~,ix] = sort(k,'descend');
    % ix = 1:numStructs;
    PlotGeneExpression(geneData(ix,iy),geneInfo(iy,:),structInfo(ix,:),true,k(ix))

case {'cerebcortex','isocortex'}
    %-------------------------------------------------------------------------------
    % Cerebral cortex enrichment: genes more strongly expressed in the cerebral cortex
    %-------------------------------------------------------------------------------
    whatTest = 'ranksum';

    % Normalize gene expression data:
    geneDataZ = BF_NormalizeMatrix(geneData,'zscore');

    switch enrichWhat
    case 'isocortex'
        isCTX = ismember(structInfo.divisionLabel,'Isocortex');
    case 'cerebCortex'
        isCTX = ismember(structInfo.divisionLabel,...
                        {'Isocortex','Olfactory Areas','Hippocampal Formation','Cortical Subplate'});
    end
    fprintf(1,'%u %s, %u non-%s\n',sum(isCTX),enrichWhat,sum(~isCTX),enrichWhat);

    gScore = zeros(numGenes,1);
    switch whatTest
    case 'ttest'
        fprintf(1,'~~TTEST P-VALUE GENE SCORES~~\n');
        for i = 1:numGenes
            [~,gScore(i)] = ttest2(geneDataZ(isCTX,i),geneDataZ(~isCTX,i),'VarType','Unequal');
        end
    case 'ranksum'
        fprintf(1,'~~RANKSUM P-VALUE GENE SCORES~~\n');
        for i = 1:numGenes
            gScore(i) = ranksum(geneDataZ(isCTX,i),geneDataZ(~isCTX,i));
        end
    end
    % Transform p-values to scores (bigger is better)
    fprintf(1,'-log10 p-values -> gene scores\n');
    gScore = -log10(gScore);

case 'genePC'
    % Genes that vary with gene expression PCs
    % Original analysis done for isocortex region filter
    geneDataNorm = BF_NormalizeMatrix(geneData,'zscore');
    fprintf(1,'Computing PCs for %ux%u matrix...',size(geneDataNorm,1),size(geneDataNorm,2));
    [pcCoeff, pcScore, ~, ~, ~] = pca(geneDataNorm,'NumComponents',2,'algorithm','als');
    fprintf(1,'\n');

    RegionScatterPlot(structInfo,pcScore(:,1),pcScore(:,2),...
                        'geneExp-PC1','geneExp-PC2','Pearson',true);

    % Score genes for correlation with the leading expression PC:
    gScore = zeros(numStructs,1);
    for i = 1:numStructs
        gScore(i) = corr(pcScore(:,1),geneData(:,i),'type','Pearson','rows','pairwise');
    end
end

if doRandomize
    warning('RANDOMIZING gene scores!!!')
    gScore = gScore(randperm(length(gScore))); % test for scores under randomization=
end

%-------------------------------------------------------------------------------
% Do the enrichment
%-------------------------------------------------------------------------------
[GOTable,geneEntrezAnnotations] = SingleEnrichment(gScore,geneInfo.entrez_id,...
                                    eParam.processFilter,eParam.sizeFilter,...
                                    eParam.numIterations);

% ANALYSIS:
numSig = sum(GOTable.pVal_corr < eParam.enrichmentSigThresh);
fprintf(1,'%u significant categories at p_corr < %.2f\n',numSig,eParam.enrichmentSigThresh);
display(GOTable(1:numSig,:));

end
