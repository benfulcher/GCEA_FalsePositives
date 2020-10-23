function [GOTable,gScore] = NodeSimpleEnrichment(params,enrichWhat,doSaveMat)
% Score each gene on some simple property then do a conventional GSR enrichment
% analysis on the results.

% ---INPUTS:
% enrichWhat = 'meanExpression'; % raw mean expression level
% enrichWhat = 'varExpression'; % raw variance of expression levels
% enrichWhat = 'cortex'; % difference in expression between cerebral cortex/others
% enrichWhat = 'genePC'; % correlation with a PC of gene expression
% enrichWhat = 'degree'; % correlation with structural connectivity degree across regions
% enrichWhat = 'excitatory'; % correlation with excitatory cell density across regions
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Set defaults:
%-------------------------------------------------------------------------------
if nargin < 1
    params = GiveMeDefaultParams('mouse');
end
if nargin < 2 || isempty(enrichWhat)
    fprintf(1,'Mean expression by default\n');
    enrichWhat = 'meanExpression'; % raw mean expression level
end
if nargin < 3
    doSaveMat = false;
end
corrType = params.e.whatCorr;
fprintf(1,'Using %s correlations (if relevant)\n',corrType);

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
numGenes = height(geneInfo);

%===============================================================================
%% -------------------------------SCORE THE GENES-------------------------------
%===============================================================================
switch enrichWhat
case {'excitatory','inhibitory','oligodendrocytes','glia','astrocytes','microglia','neurons','PV','SST','VIP'}
    %-------------------------------------------------------------------------------
    % Genes scored for correlation to excitatory cell density
    %-------------------------------------------------------------------------------
    [phenotypeVector,ia] = MatchMeCellDensity(structInfo,enrichWhat);
    structInfo = structInfo(ia,:);
    geneData = geneData(ia,:);

    gScore = zeros(numGenes,1);
    for i = 1:numGenes
        gScore(i) = corr(phenotypeVector,geneData(:,i),'type',corrType,'rows','pairwise');
    end

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
    doBinarize = true;
    [k,structInfoConn] = ComputeDegree(params,doBinarize);
    switch params.humanOrMouse
    case 'human'
        if ~all(structInfoConn.ROI_ID==structInfo.ROI_ID)
            error('Connectivity data doesn''t match gene-expression data');
        end
    case 'mouse'
        if ~all(structInfoConn.id==structInfo.id)
            error('Connectivity data doesn''t match gene-expression data');
        end
    end
    gScore = zeros(numGenes,1);
    pVals = zeros(numGenes,1);
    for i = 1:numGenes
        [gScore(i),pVals(i)] = corr(k,geneData(:,i),'type',corrType,'rows','pairwise');
    end

    % Plot:
    [~,iy] = sort(gScore,'descend');
    [~,ix] = sort(k,'descend');
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
    % How many are individually significant?:
    pCorr = mafdr(gScore,'BHFDR','true');
    fprintf(1,'%u/%u genes have corrected p < 0.05\n',sum(pCorr < 0.05),length(gScore));

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
    fprintf(1,'Scoring genes by their absolute correlation to PC1\n');
    gScore = zeros(numGenes,1);
    for i = 1:numGenes
        gScore(i) = abs(corr(pcScore(:,1),geneData(:,i),'type',corrType,'rows','pairwise'));
    end
otherwise
    error('I don''t know how to enrich for %s!',enrichWhat)
end

if doRandomize
    warning('RANDOMIZING gene scores!!!')
    gScore = gScore(randperm(length(gScore))); % test for scores under randomization=
end

%-------------------------------------------------------------------------------
% Do the enrichment
%-------------------------------------------------------------------------------
GOTable = SingleEnrichment(gScore,geneInfo.entrez_id,params.e);

%-------------------------------------------------------------------------------
% Count number of significant categories, under corrected permutation testing:
numSigPerm = sum(GOTable.pValPermCorr < params.e.sigThresh);
fprintf(1,'%u significant categories (permutation) at p_corr < %.2f\n',numSigPerm,params.e.sigThresh);
numSigZ = sum(GOTable.pValZCorr < params.e.sigThresh);
fprintf(1,'%u significant categories (gaussian-null) at p_corr < %.2f\n',numSigZ,params.e.sigThresh);
display(GOTable(1:numSigZ,:));

%-------------------------------------------------------------------------------
% Save .mat file for further analysis later
if doSaveMat
    fileSaveTo = GiveMeSimpleEnrichmentOutputFile(params,enrichWhat);
    save(fileSaveTo,'params','GOTable','gScore');
    fprintf(1,'Saved results to %s\n',fileSaveTo);
end

%-------------------------------------------------------------------------------
% Output csv for the isocortex results:
if strcmp(enrichWhat,'isocortex')
    IDLabel = GOTable.GOIDlabel;
    CategoryName = GOTable.GOName;
    ID = GOTable.GOID;
    MeanNLog10RankSumPValue = GOTable.meanScore;
    GSEA_normApprox_pFDR = GOTable.pValZCorr;
    T = table(CategoryName,IDLabel,ID,MeanNLog10RankSumPValue,GSEA_normApprox_pFDR);
    fileOut = fullfile('SupplementaryTables','MouseIsocortexRanksumGSEA.csv');
    writetable(T,fileOut,'Delimiter',',','QuoteStrings',true);
    fprintf(1,'Saved isocortex GSEA results to %s\n',fileOut);
end

end
