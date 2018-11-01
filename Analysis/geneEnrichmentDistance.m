function GOTable = geneEnrichmentDistance(params)
% Quantify for each gene the correlation between its correlated expression
% pattern and the pairwise distance between regions
% (then can try to work up a correction?)
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Inputs:
if nargin < 1
    params = GiveMeDefaultParams('mouse');
end

% Text summary of the relevant analysis settings:
textLabel = sprintf('dScores_%s_%s-%s_%s_abs-%s_conn%u',params.gcc.whatCorr,params.g.normalizationGene,...
                        params.g.normalizationRegion,params.gcc.pValOrStat,params.gcc.absType,...
                        params.gcc.onlyConnections);

%-------------------------------------------------------------------------------
% Load and process data
%-------------------------------------------------------------------------------
% Pairwise distance data:
distMat = GiveMeDistanceMatrix(params.humanOrMouse,params.c.structFilter);

% Gene expression data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
numGenes = height(geneInfo);
fprintf(1,'%u x %u gene expression matrix\n',size(geneData,1),size(geneData,2));

% Structural connectivity
if params.gcc.onlyConnections
    % (ONLY used to look specifically at connected pairs of areas: if params.gcc.onlyConnections)
    [A_bin,regionAcronyms] = GiveMeAdj(params.c.connectomeSource,...
                params.c.pThreshold,true,params.c.whatWeightMeasure,...
                params.c.whatHemispheres,params.c.structFilter);
    fprintf(1,'%u x %u structural connectivity matrix\n',size(A_bin,1),size(A_bin,2));

    % Check it aligns with the gene expression data:
    if ~all(strcmp(structInfo.acronym,regionAcronyms))
        error('Connectivity and gene expression data are not aligned');
    end
else
    A_bin = [];
end

%-------------------------------------------------------------------------------
% Filter structures:
numStructs = height(structInfo);

%-------------------------------------------------------------------------------
% Construct vector of pairwise separation distances:
dData = distMat;
if params.gcc.onlyConnections
    fprintf(1,'Only looking at pairs of regions where connections exist\n');
    dData(A_bin == 0) = 0;
else
    fprintf(1,'Looking at all pairs of regions\n');
end
% Keep upper diagonal of pairwise distances:
dData(tril(true(size(dData)),-1)) = 0;

%-------------------------------------------------------------------------------
% Visualize bulk distance trend for cge computed across all genes:
cge = corr(geneData','type','Pearson','rows','pairwise');
isUpperDiag = triu(true(size(dData)),+1);
f = figure('color','w');
plot(dData(isUpperDiag),cge(isUpperDiag),'.k');
xlabel('Distance')
ylabel('CGE (Pearson)')
title(sprintf('Bulk spatial trend: %u genes',numGenes));

%-------------------------------------------------------------------------------
% Score genes:
if params.gcc.regressDistance
    fprintf(1,'Scoring %u genes on correlation between g.gT and distance (AND regressing out distance?!)\n',numGenes);
    dRegressor = dData;
else
    fprintf(1,'Scoring %u genes on correlation between g.gT and distance (no regressor)\n',numGenes);
    dRegressor = [];
end
[geneScores,geneEntrezIDs] = GiveMeGCC(dData,geneData,geneInfo.entrez_id,params,dRegressor);
fprintf(1,'Scoring complete across %u/%u genes!\n',length(geneScores),numGenes);

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
