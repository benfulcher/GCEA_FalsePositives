function [geneData,geneInfo,structInfo] = LoadMeG(gParam,humanOrMouse)
% Load in the gene data as G from new, SDK results
% Gets the gene_energy or gene_density matrix, and the matching geneInfo table
% Along with the structInfo table.
% (could also get the full structure dataset info from the 'energy' or 'density'
% fields and match with the datasetInfo table)

%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 2
    humanOrMouse = 'mouse';
end
if nargin < 1 || isempty(gParam)
    gParam = GiveMeDefaultParams('gene',humanOrMouse);
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Get NEW DATA FROM SDK RETRIEVALS:
%-------------------------------------------------------------------------------
% try
switch humanOrMouse
case 'mouse'
    dataFile = '/Users/benfulcher/GoogleDrive/Work/CurrentProjects/CellTypesMouse/Code/Data/AllenGeneDataset_19419.mat';
    fprintf(1,'New Allen SDK-data from %s\n',dataFile);
    load(dataFile,'GeneExpData','geneInfo','structInfo');

    switch gParam.energyOrDensity
    case 'energy'
        geneData = GeneExpData.gene_energy;
    case 'density'
        geneData = GeneExpData.gene_density;
    end
case 'human'
    dataFile = 'geneSample_aparcaseg.mat'; % APARC parcellation
    % dataFile = 'geneROI_HCP.mat'; % HCP parcellation
    G = load(dataFile,'SampleGeneExpression','probeInformation');
    fprintf(1,'Loaded human expression data from %s\n',which(dataFile));
    whatROI = G.SampleGeneExpression(:,1);
    sampleExpression = G.SampleGeneExpression(:,2:end);

    % Format geneInfo into table:
    EntrezID = G.probeInformation.EntrezID;
    Symbol = G.probeInformation.GeneSymbol;
    Name = G.probeInformation.GeneName;
    DSscore = G.probeInformation.DS;
    geneInfo = table(EntrezID,Symbol,Name,DSscore);

    % Get ROI information:
    structInfo = GiveMeAPARCNames();
    % Go through ROI labels and average each:
    ROIs = sort(unique(whatROI),'ascend');
    if ~all(diff(ROIs)==1)
        error('error matching ROIs??')
    end
    numROIs = length(ROIs);
    fprintf(1,'%u ROIs\n',numROIs);
    % Match ROIs to structInfo:
    [~,ia,ib] = intersect(ROIs,structInfo.ID);
    structInfo = structInfo(ib,:);
    if length(ia)<numROIs
        error('ROIs not labeled according to APARC text file IDs...?');
    end
    fprintf(1,'Trimmed structure info to match %u structures\n',height(structInfo));
    fprintf(1,'Matching ROIs to structure info & meaning expression levels of samples in each ROI\n');
    numGenes = size(sampleExpression,2);
    geneData = zeros(numROIs,numGenes);
    for i = 1:numROIs
        geneData(i,:) = mean(sampleExpression(whatROI==ROIs(i),:),1);
    end
    fprintf(1,'Done.\n');
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
    warning('Only looking at a random set of %u genes',gParam.subsetOfGenes);
    rp = randperm(size(geneData,2));
    rp = rp(1:gParam.subsetOfGenes);
    geneData = geneData(:,rp);
    geneInfo = geneInfo(rp,:);
end

end
