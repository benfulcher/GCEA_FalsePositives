% PhenotypeNullSummary
% The idea is to understand what the nulls are likely to look like for a given phenotype

% Parameters:
whatSpecies = 'mouse';
structFilter = 'cortex';

% ---Phenotype is degree---
doBinarize = true;
params = GiveMeDefaultParams(whatSpecies,structFilter);
[k,structInfoConn] = ComputeDegree(params,doBinarize);

% Gene data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
numAreas = height(structInfo);
numGenes = height(geneInfo);

% ---NULL 1--- random-gene
nullScoreRandomGene = zeros(numGenes,1);
for i = 1:numGenes
    nullScoreRandomGene(i) = corr(k,geneData(:,i),'type','Spearman','rows','pairwise');
end
fprintf(1,'random gene: %.2g +/- %.2g\n',mean(nullScoreRandomGene),std(nullScoreRandomGene));

% ---NULL 2--- SBP-random
numNullSamples = 10000;
nullMaps = rand(numAreas,numNullSamples);
nullScoreSBPrandom = zeros(numNullSamples,1);
for i = 1:numNullSamples
    nullScoreSBPrandom(i) = corr(k,nullMaps(:,i),'type','Spearman','rows','pairwise');
end
fprintf(1,'SBPrandom: %.2g +/- %.2g\n',mean(nullScoreSBPrandom),std(nullScoreSBPrandom));

% ---NULL 3--- SBP-spatial
switch whatSpecies
case 'human'
    dataFileSurrogate = 'humanSurrogate_N20000_rho8_d02000.csv';
case 'mouse'
    dataFileSurrogate = 'mouseCortexSurrogate_N20000_rho8_d040.csv';
    % dataFileSurrogate = 'mouseSurrogate_N20000_rho8_d040.csv';;
end
nullMaps = dlmread(dataFileSurrogate,',',1,1);
nullScoreSBPspatial = zeros(numNullSamples,1);
for i = 1:numNullSamples
    nullScoreSBPspatial(i) = corr(k,nullMaps(:,i),'type','Spearman','rows','pairwise');
end
fprintf(1,'SBPspatial: %.2g +/- %.2g\n',mean(nullScoreSBPspatial),std(nullScoreSBPspatial));

%===============================================================================
f = figure('color','w');
hold('on');
histogram(nullScoreRandomGene,'normalization','pdf')
histogram(nullScoreSBPrandom,'normalization','pdf')
histogram(nullScoreSBPspatial,'normalization','pdf')
legend('randomGene','SBPrandom','SBPspatial')
title('Degree')
