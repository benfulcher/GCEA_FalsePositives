function [geneData,geneInfo,structInfo] = LoadMeG(gParam)
% Load in the gene data as G from new, SDK results
% Gets the gene_energy or gene_density matrix, and the matching geneInfo table
% Along with the structInfo table.
% (could also get the full structure dataset info from the 'energy' or 'density'
% fields and match with the datasetInfo table)

%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 1 || isempty(gParam)
    gParam = GiveMeDefaultParams('gene');
end
% normalizationSettings = {'none','none'};
% In the format: {howToNormalizeGenesAcrossRegions,howToNormalizeRegionsAcrossGenes}
% if nargin < 2
%     energyOrDensity = 'energy';
% end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Get NEW DATA FROM SDK RETRIEVALS:
%-------------------------------------------------------------------------------
% try
dataFile = '/Users/benfulcher/GoogleDrive/Work/CurrentProjects/CellTypesMouse/Code/AllenGeneDataset_19419.mat';
fprintf(1,'New Allen SDK-data from %s\n',dataFile);
load(dataFile,'GeneExpData','geneInfo','structInfo');
% catch
    % dataFile = which('AllenGeneDataset_19419.mat');
    % fprintf(1,'New Allen SDK-data from %s\n',dataFile);
    % load(dataFile,'GeneExpData','geneInfo','structInfo');
% end

switch gParam.energyOrDensity
case 'energy'
    geneData = GeneExpData.gene_energy;
case 'density'
    geneData = GeneExpData.gene_density;
end

%-------------------------------------------------------------------------------
% FURTHER PROCESSING:
%-------------------------------------------------------------------------------
geneData = BF_NormalizeMatrix(geneData,gParam.normalizationGene);
fprintf(1,'1. Normalized expression for each gene using %s\n',gParam.normalizationGene);

geneData = BF_NormalizeMatrix(geneData',gParam.normalizationRegion)';
fprintf(1,'2. Normalized expression across each brain region using %s\n',gParam.normalizationRegion);

%-------------------------------------------------------------------------------
% Subset:
%-------------------------------------------------------------------------------
if ~isempty(gParam.subsetOfGenes)
    warning('Only looking at a random set of %u genes',subsetOfGenes);
    rp = randperm(size(geneData,2));
    rp = rp(1:gParam.subsetOfGenes);
    geneData = geneData(:,rp);
    geneInfo = geneInfo(rp,:);
end

end
