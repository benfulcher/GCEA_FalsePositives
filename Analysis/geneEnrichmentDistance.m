function geneEnrichmentDistance()
% Quantify for each gene the correlation between its coexpression and the
% pairwise distance between regions
%-------------------------------------------------------------------------------
% Then can try to work up a correction
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Parameters:
%-------------------------------------------------------------------------------
whatCorr = 'Spearman';
energyOrDensity = 'energy'; % what gene expression data to use
pValOrStat = 'stat'; % 'pval','stat'

%-------------------------------------------------------------------------------
% Load in and process data:
%-------------------------------------------------------------------------------
% Connectivity data:
C = load('Mouse_Connectivity_Data.mat');

% Gene expression data:
G = LoadMeG();
GData = G.GeneExpData.(energyOrDensity);
% z-score across genes for each brain region:
zGData = BF_NormalizeMatrix(GData','zscore')';

% Pairwise euclidean distance data for control:
d = C.Dist_Matrix{1,1}/1000;
d_upper = d(triu(true(size(d)),1));

%-------------------------------------------------------------------------------
% Now we want to score genes:
%-------------------------------------------------------------------------------
rsfMRI_upper = rsfMRI(triu(true(size(rsfMRI)),1));
for i = 1:G.numGenes
    % Compute this gene's correlation pattern:
    g = zGData(:,i);
    % gz = (g-nanmean(g))/nanstd(g);
    ggBlock = g*g';
    ggBlock_upper = ggBlock(triu(true(size(ggBlock)),1));

    % Only compute correlation for gene score when the proportion of missing
    % data is under 20%:
    if mean(isnan(ggBlock_upper)) < 0.2
        isGood = ~isnan(ggBlock_upper);
        [rho,pVal] = corr(ggBlock_upper(isGood),d_upper(isGood),'type',whatCorr);
        % [rho,pVal] = corr(ggBlock_upper(isGood),rsfMRI_upper(isGood),'type',whatCorr);
        switch pValOrStat
        case 'pVal'
            geneScores(i) = pVal;
        case 'stat'
            geneScores(i) = abs(rho);
        end
    else
        geneScores(i) = NaN;
    end
end


%-------------------------------------------------------------------------------
% Plot the genes with the highest scores, or a range of scores
%-------------------------------------------------------------------------------
% [~,ix] = sort(geneScores,'ascend');
% f = figure('color','w');
% plot()

%-------------------------------------------------------------------------------
% Now we can run enrichment on them:
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
% Write them out to ermine J:
writeErmineJFile(sprintf('corr_d_%s_%s_%s',whatCorr,whatGeneData,pValOrStat),...
                geneScoresWrite,geneEntrezWrite,sprintf('%s_%s_%s',whatCorr,whatGeneData,pValOrStat));

end
