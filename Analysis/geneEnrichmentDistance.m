function GOTable = geneEnrichmentDistance(structFilter,whatSpecies)
% Quantify for each gene the correlation between its coexpression and the
% pairwise distance between regions
%-------------------------------------------------------------------------------
% Then can try to work up a correction
%-------------------------------------------------------------------------------

if nargin < 1 || isempty(structFilter)
    structFilter = 'isocortex'; % 'cortex', 'all'
end
if nargin < 2
    whatSpecies = 'mouse'; % 'mouse', 'human'
end

%-------------------------------------------------------------------------------
% Gene scoring parameters:
%-------------------------------------------------------------------------------
whatCorr = 'Spearman'; % 'Pearson', 'Spearman'
pValOrStat = 'stat'; % 'pval','stat'
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
absType = 'neg';
onlyConnections = false; % only look where there are connections
% regressDistance = true; % whether to regress distance

%-------------------------------------------------------------------------------
% Get default parameter sets:
params = GiveMeDefaultParams(whatSpecies);

%-------------------------------------------------------------------------------
% Load in and process data:
%-------------------------------------------------------------------------------
% Distance data:
distMat = GiveMeDistanceMatrix(whatSpecies);

% Binary connectome:
[A_bin,regionAcronyms,adjPVals] = GiveMeAdj(params.c.connectomeSource,...
            params.c.pThreshold,true,params.c.whatWeightMeasure,...
            params.c.whatHemispheres,params.c.structFilter);

%-------------------------------------------------------------------------------
% Gene expression data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
numGenes = height(geneInfo);

%-------------------------------------------------------------------------------
% Filter to subregions
%-------------------------------------------------------------------------------
% Filter structures:
[A_bin,geneData,structInfo,keepStruct] = filterStructures(structFilter,structInfo,A_bin,geneData);
distMat = distMat(keepStruct,keepStruct);
entrezIDs = geneInfo.entrez_id;
numStructs = height(structInfo);

%-------------------------------------------------------------------------------
% Make vector of distances
%-------------------------------------------------------------------------------
if onlyConnections
    pThreshold = 0.05;
    dData = distMat;
    dData(A_bin == 0) = 0;
    % Only symmetric?:
    dData(tril(true(size(dData)),-1)) = 0;
else
    % Convert to vector of upper diagonal
    dData = distMat;
    dData(tril(true(size(dData)),-1)) = 0;
end

%-------------------------------------------------------------------------------
% Score genes:
%-------------------------------------------------------------------------------
% if regressDistance
%     fprintf(1,'Scoring %u genes on coexpression with distance (AND regressing out distance?!)\n',numGenes);
%     [geneScores,geneEntrezIDs] = GiveMeGCC(dData,geneData,entrezIDs,whatCorr,...
%                                 dData,absType,thresholdGoodGene,pValOrStat);
% else
fprintf(1,'Scoring %u genes on coexpression with distance\n',numGenes);
[geneScores,geneEntrezIDs] = GiveMeGCC(dData,geneData,entrezIDs,whatCorr,...
                            [],absType,thresholdGoodGene,pValOrStat);
% end
fprintf(1,'Scoring done. Enrichment time\n');

%-------------------------------------------------------------------------------
% Do the enrichment:
%-------------------------------------------------------------------------------
GOTable = SingleEnrichment(geneScores,geneEntrezIDs,...
                                    params.e.whatSource,...
                                    params.e.processFilter,params.e.sizeFilter,...
                                    params.e.numIterations);

% ANALYSIS:
numSig = sum(GOTable.pValCorr < params.e.enrichmentSigThresh);
fprintf(1,'%u significant categories at p_corr < %.2f\n',numSig,params.e.enrichmentSigThresh);
display(GOTable(1:numSig,:));

%-------------------------------------------------------------------------------
textLabel = sprintf('dScores_%s_%s-%s_%s_abs-%s_conn%u',whatCorr,params.g.normalizationGene,...
                        params.g.normalizationRegion,pValOrStat,absType,onlyConnections);

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
% geneEntrez = geneEntrezIDs;
% geneDistanceScores = geneScores;
% fileNameMat = [textLabel,'.mat'];
% save(fileNameMat,'geneEntrez','geneDistanceScores');
% fprintf(1,'Saved info to %s\n',fileNameMat);

%-------------------------------------------------------------------------------
% Do enrichment using ermine J:
%-------------------------------------------------------------------------------
% numIterations = 20000;
% fileNameWrite = writeErmineJFile('tmp',-(geneScores),geneEntrezIDs,'distance');
% ermineJResults = RunErmineJ(fileNameWrite,numIterations);

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

end
