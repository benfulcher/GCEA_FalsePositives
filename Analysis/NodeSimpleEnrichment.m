function [GOTable,gScore] = NodeSimpleEnrichment(params,enrichWhat,structFilter,corrType)
% Score each gene on some simple property

% ---INPUTS:
% enrichWhat = 'meanExpression'; % raw mean expression level
% enrichWhat = 'varExpression'; % raw variance of expression levels
% enrichWhat = 'cortex'; % difference in expression between cerebral cortex/others
% enrichWhat = 'genePC'; % correlation with a PC of gene expression
% enrichWhat = 'degree'; % correlation with structural connectivity degree across regions
%
% structFilter = 'cortex';
%
% corrType = 'Pearson';
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Name an analysis, and this script loads the data and does the enrichment
%-------------------------------------------------------------------------------
if nargin < 1
    params = GiveMeDefaultParams('mouse');
end
if nargin < 2 || isempty(enrichWhat)
    fprintf(1,'Mean expression by default\n');
    enrichWhat = 'meanExpression'; % raw mean expression level
end
if nargin < 3 || isempty(structFilter)
    fprintf(1,'No filter: all brain regions included\n');
    structFilter = 'all'; % 'all', 'isocortex'
end
if nargin < 4 || isempty(corrType)
    corrType = 'Pearson';
    fprintf(1,'Pearson correlations by default (if relevant)\n');
end

doRandomize = false;

%-------------------------------------------------------------------------------
% Warn if normalization applied:
if ismember(enrichWhat,{'meanExpression','varExpression'})
    if ~strcmp(params.g.normalizationGene,'none')
        warning('Gene expression normalization across genes?? (Consider -> ''none'') :-O')
    end
    if ~strcmp(params.g.normalizationRegion,'none')
        warning('Gene expression normalization across regions?? (Consider ->''none'') :-O')
    end
end

%-------------------------------------------------------------------------------
% Gene data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
% Filter structures:
[~,geneData,structInfo] = filterStructures(structFilter,structInfo,[],geneData);
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
    % Degree enrichment: genes scored for the correlation of their expression to
    % the number of connections each region makes to other regions
    %-------------------------------------------------------------------------------
    % Binary connectome data:
    k = ComputeDegree(params.humanOrMouse,true);
    gScore = zeros(numGenes,1);
    for i = 1:numGenes
        gScore(i) = corr(k,geneData(:,i),'type',corrType,'rows','pairwise');
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
    fprintf(1,'Using a ranksum test\n');

    % Normalize gene expression data:
    geneDataZ = BF_NormalizeMatrix(geneData,'zscore');
    warning('Ensuring z-score normalized gene expression data');

    switch params.humanOrMouse
    case 'mouse'
        switch enrichWhat
        case 'isocortex'
            isCTX = ismember(structInfo.divisionLabel,'Isocortex');
        case 'cerebCortex'
            isCTX = ismember(structInfo.divisionLabel,...
                            {'Isocortex','Olfactory Areas','Hippocampal Formation','Cortical Subplate'});
        end
    case 'human'
        isCTX = structInfo.isCortex;
    end
    fprintf(1,'%u %s, %u non-%s\n',sum(isCTX),enrichWhat,sum(~isCTX),enrichWhat);
    if ~any(isCTX) || ~any(~isCTX)
        error('Not enough data to distinguish cortical from non-cortical areas');
    end

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
    fprintf(1,'Computing PCs for %ux%u matrix using the als algorithm...',...
                                size(geneDataNorm,1),size(geneDataNorm,2));
    [pcCoeff, pcScore, ~, ~, ~] = pca(geneDataNorm,'NumComponents',2,'algorithm','als');
    fprintf(1,'\n');

    if strcmp(params.humanOrMouse,'mouse')
        RegionScatterPlot(structInfo,pcScore(:,1),pcScore(:,2),...
                        'geneExp-PC1','geneExp-PC2',corrType,true);
    end

    % Score genes for correlation with the leading expression PC:
    fprintf(1,'Scoring genes by their correlation to PC1\n');
    gScore = zeros(numGenes,1);
    for i = 1:numGenes
        gScore(i) = corr(pcScore(:,1),geneData(:,i),'type',corrType,'rows','pairwise');
    end
end

if doRandomize
    warning('RANDOMIZING gene scores!!!')
    gScore = gScore(randperm(length(gScore))); % test for scores under randomization=
end

%-------------------------------------------------------------------------------
% Do the enrichment
%-------------------------------------------------------------------------------
GOTable = SingleEnrichment(gScore,geneInfo.entrez_id,params.e);

% ANALYSIS:
numSig = sum(GOTable.pValCorr < params.e.sigThresh);
fprintf(1,'%u significant categories at p_corr < %.2f\n',numSig,params.e.sigThresh);
display(GOTable(1:numSig,:));

end
