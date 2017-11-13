function GOTable = geneEnrichmentDistance(structFilter,whatSpecies,params,GCCparams)
% Quantify for each gene the correlation between its coexpression and the
% pairwise distance between regions
%-------------------------------------------------------------------------------
% Then can try to work up a correction
%-------------------------------------------------------------------------------

if nargin < 1 || isempty(structFilter)
    structFilter = 'isocortex'; % 'cortex', 'all'
end
if nargin < 2 || isempty(whatSpecies)
    whatSpecies = 'mouse'; % 'mouse', 'human'
end
if nargin < 3
    params = GiveMeDefaultParams(whatSpecies);
end
if nargin < 4
    GCCparams = struct();
    GCCparams.whatCorr = 'Spearman'; % 'Pearson', 'Spearman'
    GCCparams.pValOrStat = 'stat'; % 'pval','stat'
    GCCparams.thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
    GCCparams.absType = 'neg';
    GCCparams.onlyConnections = false; % only look where there are structural connections
    GCCparams.regressDistance = false; % whether to regress distance
end

%-------------------------------------------------------------------------------
% Gene scoring parameters:
%-------------------------------------------------------------------------------

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
if GCCparams.onlyConnections
    fprintf(1,'Only looking at pairs of regions where connections exist\n');
    pThreshold = 0.05;
    dData = distMat;
    dData(A_bin == 0) = 0;
    % Only symmetric?:
    dData(tril(true(size(dData)),-1)) = 0;
else
    fprintf(1,'Looking at all pairs of regions\n');
    % Convert to vector of upper diagonal
    dData = distMat;
    dData(tril(true(size(dData)),-1)) = 0;
end

%-------------------------------------------------------------------------------
% Score genes:
%-------------------------------------------------------------------------------
if GCCparams.regressDistance
    fprintf(1,'Scoring %u genes on coexpression with distance (AND regressing out distance?!)\n',numGenes);
    dRegressor = dData;
else
    fprintf(1,'Scoring %u genes on coexpression with distance (no regressor)\n',numGenes);
    dRegressor = [];
end
[geneScores,geneEntrezIDs] = GiveMeGCC(dData,geneData,entrezIDs,GCCparams.whatCorr,...
                            dRegressor,GCCparams.absType,GCCparams.thresholdGoodGene,...
                            GCCparams.pValOrStat);
fprintf(1,'Gene scoring done across %u/%u genes! Enrichment time!\n',length(geneScores),numGenes);

%-------------------------------------------------------------------------------
% Do the enrichment:
%-------------------------------------------------------------------------------
GOTable = SingleEnrichment(geneScores,geneEntrezIDs,...
                                params.e.whatSource,params.e.processFilter,...
                                params.e.sizeFilter,params.e.numIterations);

% ANALYSIS:
numSig = sum(GOTable.pValCorr < params.e.enrichmentSigThresh);
fprintf(1,'%u significant categories at p_corr < %.2f\n',numSig,params.e.enrichmentSigThresh);
display(GOTable(1:numSig,:));

%-------------------------------------------------------------------------------
textLabel = sprintf('dScores_%s_%s-%s_%s_abs-%s_conn%u',GCCparams.whatCorr,params.g.normalizationGene,...
                        params.g.normalizationRegion,GCCparams.pValOrStat,GCCparams.absType,...
                        GCCparams.onlyConnections);

%-------------------------------------------------------------------------------
% Look at the distribution
%-------------------------------------------------------------------------------
f = figure('color','w');
histogram(geneScores)
xlabel(sprintf('%s correlation between distance and coexpression',GCCparams.whatCorr))
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
