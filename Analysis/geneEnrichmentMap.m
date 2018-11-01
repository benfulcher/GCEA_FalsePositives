function GOTable = geneEnrichmentMap(params,myMap)
% Quantify for each gene the correlation between its correlated expression
% pattern and the pairwise distance between regions
% (then can try to work up a correction?)
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Inputs:
if nargin < 1
    params = GiveMeDefaultParams('mouse');
end

%-------------------------------------------------------------------------------
% Load and process data
%-------------------------------------------------------------------------------
% Gene expression data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
fprintf(1,'%u x %u gene expression matrix\n',size(geneData,1),size(geneData,2));
numGenes = height(geneInfo);
numStructs = height(structInfo);

%-------------------------------------------------------------------------------
% Compute gene scores as correlations to the given spatial map


%-------------------------------------------------------------------------------
% Do the enrichment:
fprintf(1,'Enrichment time!\n');
GOTable = SingleEnrichment(geneScores,geneEntrezIDs,params.e);

%-------------------------------------------------------------------------------
% Look at the distribution of scores across genes and within categories
f = figure('color','w');
hold('on')
h1 = histogram(geneScores);
h2 = histogram(GOTable.meanScore);
xlabel(sprintf('%s correlation between CGE and separation distance',params.gcc.whatCorr))
legend([h1,h2],{sprintf('%u gene scores',length(geneScores)),sprintf('%u category scores',height(GOTable))})
title(textLabel,'interpreter','none')

%-------------------------------------------------------------------------------
% List GO categories with significant p-values:
numSig = sum(GOTable.pValCorr < params.e.sigThresh);
fprintf(1,'%u GO categories have p_corr < %.2f\n',numSig,params.e.sigThresh);
display(GOTable(1:numSig,:));

%-------------------------------------------------------------------------------
% Save result to .mat file
%-------------------------------------------------------------------------------
% geneEntrez = geneEntrezIDs;
% geneDistanceScores = geneScores;
% fileNameMat = [textLabel,'.mat'];
% save(fileNameMat,'geneEntrez','geneDistanceScores');
% fprintf(1,'Saved info to %s\n',fileNameMat);

end
