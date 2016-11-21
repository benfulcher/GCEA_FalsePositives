% function geneEnrichmentDistance()
% Quantify for each gene the correlation between its coexpression and the
% pairwise distance between regions
%-------------------------------------------------------------------------------
% Then can try to work up a correction
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Parameters:
%-------------------------------------------------------------------------------
whatCorr = 'Spearman'; % 'Pearson', 'Spearman'
energyOrDensity = 'energy'; % what gene expression data to use
pValOrStat = 'stat'; % 'pval','stat'
normalizeHow = 'raw'; % 'mixedSigmoid', 'zscore'
thresholdGoodGene = 0.2; % threshold of valid coexpression values at which a gene is kept

%-------------------------------------------------------------------------------
% Load in and process data:
%-------------------------------------------------------------------------------
% Connectivity data:
C = load('Mouse_Connectivity_Data.mat');

%-------------------------------------------------------------------------------
% Gene expression data:
G = LoadMeG();
GData.raw = G.GeneExpData.(energyOrDensity);
numGenes = size(GData.raw,2);
% Normalize expression levels across brain regions for each gene:
% GData.z = BF_NormalizeMatrix(GData.raw,normalizeHow);
% z-score across genes for each brain region:
% GData.zz = BF_NormalizeMatrix(GData.z','zscore')';
% fprintf(1,'Gene data normalized.\n');

%-------------------------------------------------------------------------------
% Pairwise Euclidean distance data:
d = C.Dist_Matrix{1,1}/1000;
% Convert to vector of upper diagonal
d_upper = d(triu(true(size(d)),1));

%-------------------------------------------------------------------------------
% Score genes:
%-------------------------------------------------------------------------------
fprintf(1,'Scoring %u genes on coexpression with distance\n',numGenes);
for i = 1:numGenes
    % This gene's correlation pattern across regions:
    g = GData.raw(:,i);
    % Correlation with itself (assuming normal distribution of expression):
    ggBlock = g*g';
    % Convert to vector of upper diagonal
    ggBlock_upper = ggBlock(triu(true(size(ggBlock)),1));

    % Only compute correlation for gene score when the proportion of missing
    % data is under 20%:
    if mean(isnan(ggBlock_upper)) < thresholdGoodGene
        [rho,pVal] = corr(ggBlock_upper,d_upper,'rows','pairwise','type',whatCorr);
        switch pValOrStat
        case 'pVal'
            geneScores(i) = pVal;
        case 'stat'
            geneScores(i) = rho;
        end
    else
        geneScores(i) = NaN;
    end
end

%-------------------------------------------------------------------------------
% Filter to genes with valid scores:
%-------------------------------------------------------------------------------
geneEntrez = [G.GeneStruct.gene_entrez_id];
% Filtering:
if any(~isfinite(geneScores))
    keepGenes = isfinite(geneScores);
    fprintf(1,'Removing %u genes -> now %u\n',sum(~keepGenes),sum(keepGenes));
    geneScoresWrite = geneScores(keepGenes);
    geneEntrezWrite = geneEntrez(keepGenes);
else
    geneScoresWrite = geneScores;
    geneEntrezWrite = geneEntrez;
end

%-------------------------------------------------------------------------------
% Look at the distribution
%-------------------------------------------------------------------------------
f = figure('color','w');
histogram(geneScoresWrite)
xlabel(sprintf('%s correlation between distance and coexpression',whatCorr))

%-------------------------------------------------------------------------------
% Save result to .mat file
%-------------------------------------------------------------------------------
fileNameMat = sprintf('dScores_%s.mat',whatCorr);
geneEntrez = geneEntrezWrite;
geneDistanceScores = geneScoresWrite;
save(fileNameMat,'geneEntrez','geneDistanceScores');

%-------------------------------------------------------------------------------
% Write them out to ermine J:
%-------------------------------------------------------------------------------
fileNameBase = sprintf('corr_d_%s_%s_%s',whatCorr,normalizeHow,pValOrStat);
doWhat = {'pos','neg','abs'};
for i = 1:length(doWhat)
    fileName = sprintf('%s_%s',fileNameBase,doWhat{i});
    switch doWhat{i}
    case 'pos'
        geneScoresNow = geneScoresWrite;
    case 'neg'
        geneScoresNow = -geneScoresWrite;
    case 'abs'
        geneScoresNow = abs(geneScoresWrite);
    end

    writeErmineJFile(fileName,geneScoresNow,geneEntrezWrite,...
                    sprintf('%s_%s_%s',whatCorr,normalizeHow,pValOrStat));
end

% end
