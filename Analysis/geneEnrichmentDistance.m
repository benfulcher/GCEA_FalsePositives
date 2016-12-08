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
normalizationSettings = {'none','none'}; % normalization by gene/region
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
doAbs = false;
onlyConnections = false; % only look where there are connections

%-------------------------------------------------------------------------------
% Load in and process data:
%-------------------------------------------------------------------------------
% Connectivity data:
C = load('Mouse_Connectivity_Data.mat');
% Pairwise Euclidean distance data:
d = C.Dist_Matrix{1,1}/1000;

% Make vector of distances
if onlyConnections
    pThreshold = 0.05;
    A_bin = GiveMeAdj(C,'binary','ipsi',0,pThreshold);
    dVector = d(A_bin > 0);
else
    % Convert to vector of upper diagonal
    dVector = d(triu(true(size(d)),1));
end

%-------------------------------------------------------------------------------
% Gene expression data:
[geneData,geneInfo,structInfo] = LoadMeG(normalizationSettings,energyOrDensity);
numGenes = size(geneData,2);

%-------------------------------------------------------------------------------
% Score genes:
%-------------------------------------------------------------------------------
fprintf(1,'Scoring %u genes on coexpression with distance\n',numGenes);

[geneScores,geneEntrezIDs] = GiveMeGCC(dVector,geneData,geneInfo.entrez_id,whatCorr,...
                                false,doAbs,thresholdGoodGene,pValOrStat);

textLabel = sprintf('dScores_%s_%s-%s_%s_abs%u_conn%u',whatCorr,normalizationSettings{1},...
                        normalizationSettings{2},pValOrStat,doAbs,onlyConnections);
%-------------------------------------------------------------------------------
% Look at the distribution
%-------------------------------------------------------------------------------
f = figure('color','w');
histogram(geneScores)
xlabel(sprintf('%s correlation between distance and coexpression',whatCorr))
title(textLabel,'interpreter','none')

%-------------------------------------------------------------------------------
% Save result to .mat file
%-------------------------------------------------------------------------------
geneEntrez = geneEntrezIDs;
geneDistanceScores = geneScores;
save([textLabel,'.mat'],'geneEntrez','geneDistanceScores');
fprintf(1,'Saved info to %s\n',fileNameMat);

%-------------------------------------------------------------------------------
% Do enrichment using ermine J:
%-------------------------------------------------------------------------------
numIterations = 20000;
fileNameWrite = writeErmineJFile('tmp',-(geneScores),geneEntrezIDs,'distance');
% fileNameWrite = writeErmineJFile('tmp',-geneDistanceScores,geneEntrez,'distance');
ermineJResults = RunErmineJ(fileNameWrite,numIterations);

% fileName = sprintf('corr_d_%s_%s_%s',whatCorr,normalizeHow,pValOrStat);
%
% doWhat = {'pos','neg','abs'};
% for i = 1:length(doWhat)
%     fileName = sprintf('%s_%s',fileNameBase,doWhat{i});
%     switch doWhat{i}
%     case 'pos'
%         geneScoresNow = geneScores;
%     case 'neg'
%         geneScoresNow = -geneScores;
%     case 'abs'
%         geneScoresNow = abs(geneScores);
%     end
%
%     writeErmineJFile(fileName,geneScoresNow,geneEntrez,...
%                     sprintf('%s_%s_%s',whatCorr,normalizeHow,pValOrStat));
% end

% end
