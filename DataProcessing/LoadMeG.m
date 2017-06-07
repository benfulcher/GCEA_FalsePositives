function [GeneExpData,geneInfo,structInfo] = LoadMeG(normalizationSettings,energyOrDensity)
% Load in the gene data as G from new, SDK results
% Gets the gene_energy or gene_density matrix, and the matching geneInfo table
% Along with the structInfo table.
% (could also get the full structure dataset info from the 'energy' or 'density'
% fields and match with the datasetInfo table)

%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 1 || isempty(normalizationSettings)
    normalizationSettings = {'none','none'};
    % In the format: {howToNormalizeGenesAcrossRegions,howToNormalizeRegionsAcrossGenes}
end
if nargin < 2
    energyOrDensity = 'energy';
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Get NEW DATA FROM SDK RETRIEVALS:
%-------------------------------------------------------------------------------
try
    dataFile = '/Users/benfulcher/GoogleDrive/Work/CurrentProjects/CellTypesMouse/Code/AllenGeneDataset_19419.mat';
    load(dataFile,'GeneExpData','geneInfo','structInfo');
    fprintf(1,'New Allen SDK-data from %s\n',dataFile);
catch
    dataFile = which('AllenGeneDataset_19419.mat');
    fprintf(1,'New Allen SDK-data from %s\n',dataFile);
    load(dataFile,'GeneExpData','geneInfo','structInfo');
end

switch energyOrDensity
case 'energy'
    GeneExpData = GeneExpData.gene_energy;
case 'density'
    GeneExpData = GeneExpData.gene_density;
end

%-------------------------------------------------------------------------------
% FURTHER PROCESSING:
%-------------------------------------------------------------------------------
GeneExpData = BF_NormalizeMatrix(GeneExpData,normalizationSettings{1});
fprintf(1,'1. Normalized expression for each gene using %s\n',normalizationSettings{1});

GeneExpData = BF_NormalizeMatrix(GeneExpData',normalizationSettings{2})';
fprintf(1,'2. Normalized expression across each brain region using %s\n',normalizationSettings{2});

end
