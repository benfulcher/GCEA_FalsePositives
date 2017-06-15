%-------------------------------------------------------------------------------
% Name an analysis, and this script loads the data and does the enrichment
%-------------------------------------------------------------------------------
enrichWhat = 'meanExpression';

%-------------------------------------------------------------------------------
% Connectome parameters
%-------------------------------------------------------------------------------
structureFilter = 'cortex';

% Binary connectome data:
cParam = GiveMeDefaultParams('conn');
[A_bin,regionAcronyms,adjPVals] = GiveMeAdj(cParam.connectomeSource,cParam.pThreshold,true,...
                                    cParam.whatWeightMeasure,cParam.whatHemispheres,cParam.structFilter);
% Gene data:
gParam = GiveMeDefaultParams('gene');
[geneData,geneInfo,structInfo] = LoadMeG({gParam.normalizationGene,gParam.normalizationRegion},gParam.energyOrDensity);

% Filter structures:
[A_bin,geneData,structInfo] = filterStructures(structFilter,structInfo,A_bin,geneData)
numStructs = height(structInfo);

%===============================================================================
%% -------------------------------SCORE THE GENES-------------------------------
%===============================================================================
switch enrichWhat
case 'meanExpression'
    %-------------------------------------------------------------------------------
    % Mean expression enrichment: genes scored for mean expression across the brain
    %-------------------------------------------------------------------------------
    gScore = nanmean(geneData,1);
    % gScore = gScore(randperm(length(gScore))); % test for scores under randomization=

case 'varExpression'
    %-------------------------------------------------------------------------------
    % Expression variance enrichment: genes scored for variance of expression across brain
    %-------------------------------------------------------------------------------
    gScore = nanstd(geneData,1);
    % gScore = gScore(randperm(length(gScore)));

case 'degree'
    %-------------------------------------------------------------------------------
    % Degree enrichment: genes scored for number of connections they make to other regions
    %-------------------------------------------------------------------------------
    k = sum(A_bin,1)' + sum(A_bin,2);
    gScore = zeros(height(geneInfo),1);
    for i = 1:height(geneInfo)
        gScore(i) = corr(k,geneData(:,i),'type','Pearson','rows','pairwise');
    end
    % Plot:
    [~,iy] = sort(gScore,'descend');
    [~,ix] = sort(k,'descend');
    % ix = 1:numStructs;
    PlotGeneExpression(geneData(ix,iy),geneInfo(iy,:),structInfo(ix,:),true,k(ix))

case 'cortex'
    %-------------------------------------------------------------------------------
    % Cerebral cortex enrichment: genes more strongly expressed in the cerebral cortex
    %-------------------------------------------------------------------------------
    % Normalize gene expression data:
    geneDataZ = BF_NormalizeMatrix(geneData,'zscore');

    gScore = zeros(height(geneInfo),1);
    isCTX = ismember(structInfo.divisionLabel,...
                    {'Isocortex','Olfactory Areas','Hippocampal Formation','Cortical Subplate'});
    for i = 1:height(geneInfo)
        gScore(i) = ttest2(geneDataZ(isCTX,i),geneDataZ(~isCTX,i),'VarType','Unequal');
    end

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
    gScore = zeros(height(geneInfo),1);
    for i = 1:height(geneInfo)
        gScore(i) = corr(pcScore(:,1),geneData(:,i),'type','Pearson','rows','pairwise');
    end

end

%-------------------------------------------------------------------------------
% Do the enrichment
%-------------------------------------------------------------------------------
[GOTable,geneEntrezAnnotations] = SingleEnrichment(gScore,geneInfo.entrez_id,'biological_process',[5,200],20000);

% ANALYSIS:
numSig = sum(GOTable.pVal_corr < 0.05);
fprintf(1,'%u significant categories at p_corr < 0.05\n',numSig);
display(GOTable(1:numSig,:));
