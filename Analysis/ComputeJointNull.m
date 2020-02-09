function nullDistribution = ComputeJointNull(params,categorySize,numSamples)
% Joint null for random collection of genes to a given ensemble
%-------------------------------------------------------------------------------
if nargin < 2
    categorySize = 40;
end
if nargin < 3
    numSamples = 1000;
end
%-------------------------------------------------------------------------------
% Load gene-expression data:
[geneDataReal,geneInfoReal,structInfoReal] = LoadMeG(params.g);
numGenes = height(geneInfoReal);
numAreas = height(structInfoReal);

%-------------------------------------------------------------------------------
% Load ensemble data as 'nullMaps'
switch params.e.whatEnsemble
case 'randomMap'
    nullMaps = rand(numAreas,numSamples);
case 'customEnsemble'
    load(params.e.dataFileSurrogate,'nullMaps');
end

%-------------------------------------------------------------------------------
% Follow similar to PermuteForNullDistributions:
nullDistribution = zeros(1,numSamples);
fprintf(1,'Computing a joint null of ensemble--categories containing %u random genes\n',...
                    categorySize);
for j = 1:numSamples
    % Take a random phenotype:
    myPhenotype = nullMaps(:,j);
    % Take a random selection of genes:
    rp = randperm(numGenes);
    geneScores = zeros(categorySize,1);
    for i = 1:categorySize
        geneScores(i) = corr(myPhenotype,geneDataReal(:,rp(i)),...
                                        'type','Spearman','rows','pairwise');
    end
    nullDistribution(j) = nanmean(geneScores);
end

end
