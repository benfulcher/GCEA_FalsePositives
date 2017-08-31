function geneEnrichmentDistance(structFilter)
% Quantify for each gene the correlation between its coexpression and the
% pairwise distance between regions
%-------------------------------------------------------------------------------
% Then can try to work up a correction
%-------------------------------------------------------------------------------

if nargin < 1 || isempty(structFilter)
    structFilter = 'cortex'; % 'cortex', 'all'
end

%-------------------------------------------------------------------------------
% Gene scoring parameters:
%-------------------------------------------------------------------------------
whatCorr = 'Spearman'; % 'Pearson', 'Spearman'
pValOrStat = 'stat'; % 'pval','stat'
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept
absType = 'neg';
onlyConnections = false; % only look where there are connections

%-------------------------------------------------------------------------------
% Get default parameter sets:
cParam = GiveMeDefaultParams('conn');
eParam = GiveMeDefaultParams('enrichment');
gParam = GiveMeDefaultParams('gene');
gParam.normalizationGene = 'none'; % 'none', 'mixedSigmoid'
gParam.normalizationRegion = 'none'; % 'none', 'zscore'

%-------------------------------------------------------------------------------
% Load in and process data:
%-------------------------------------------------------------------------------
% Distance data:
C = load('Mouse_Connectivity_Data.mat','Dist_Matrix');
% Pairwise ipsilateral Euclidean distance data:
d = C.Dist_Matrix{1,1}/1000;
% Binary connectome:
[A_bin,regionAcronyms,adjPVals] = GiveMeAdj(cParam.connectomeSource,cParam.pThreshold,true,...
                                cParam.whatWeightMeasure,cParam.whatHemispheres,cParam.structFilter);

%-------------------------------------------------------------------------------
% Gene expression data:
[geneData,geneInfo,structInfo] = LoadMeG(gParam);

%-------------------------------------------------------------------------------
% Filter to subregions
%-------------------------------------------------------------------------------
% Filter structures:
[A_bin,geneData,structInfo] = filterStructures(structFilter,structInfo,A_bin,geneData);
if strcmp(structFilter,'cortex')
    keepStruct = strcmp(structInfo.divisionLabel,'Isocortex');
    geneData = geneData(keepStruct,:);
    structInfo = structInfo(keepStruct,:);
    A_bin = A_bin(keepStruct,keepStruct);
    d = d(keepStruct,keepStruct);
end
entrezIDs = geneInfo.entrez_id;
numGenes = height(geneInfo);
numStructs = height(structInfo);

%-------------------------------------------------------------------------------
% Make vector of distances
%-------------------------------------------------------------------------------
if onlyConnections
    pThreshold = 0.05;
    dData = d;
    dData(A_bin == 0) = 0;
    % Only symmetric?:
    dData(tril(true(size(dData)),-1)) = 0;
else
    % Convert to vector of upper diagonal
    dData = d;
    dData(tril(true(size(d)),-1)) = 0;
end

%-------------------------------------------------------------------------------
% Score genes:
%-------------------------------------------------------------------------------
fprintf(1,'Scoring %u genes on coexpression with distance\n',numGenes);

[geneScores,geneEntrezIDs] = GiveMeGCC(dData,geneData,entrezIDs,whatCorr,...
                                [],absType,thresholdGoodGene,pValOrStat);

%-------------------------------------------------------------------------------
% Do the enrichment:
%-------------------------------------------------------------------------------
[GOTable,geneEntrezAnnotations] = SingleEnrichment(geneScores,geneEntrezIDs,...
                                    eParam.processFilter,eParam.sizeFilter,...
                                    eParam.numIterations);

% ANALYSIS:
numSig = sum(GOTable.pVal_corr < eParam.enrichmentSigThresh);
fprintf(1,'%u significant categories at p_corr < %.2f\n',numSig,eParam.enrichmentSigThresh);
display(GOTable(1:numSig,:));

%-------------------------------------------------------------------------------
textLabel = sprintf('dScores_%s_%s-%s_%s_abs-%s_conn%u',whatCorr,gParam.normalizationGene,...
                        gParam.normalizationRegion,pValOrStat,absType,onlyConnections);

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
