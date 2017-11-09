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
    gParam = GiveMeDefaultParams('gene',humanOrMouse);
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Get NEW DATA FROM SDK RETRIEVALS:
%-------------------------------------------------------------------------------
% try
switch gParam.humanOrMouse
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
    % Piece together filename from parameters:
    switch gParam.whatParcellation
    case 'APARC'
        dataFileBase = 'aparcaseg_'; % APARC parcellation
    case 'HCP'
        dataFileBase = '360parcellationLcortex_'; % HCP parcellation
    end
    switch gParam.probeSelection
    case 'mean'
        dataFileProbe = '_ProbeMean';
    case 'variance'
        dataFileProbe = '_ProbeVariance';
    end
    switch gParam.normalizationInternal
    case 'robustSigmoid'
        dataFileNorm = '_RobustSigmoid';
    case 'none'
        dataFileNorm = '';
    end
    % Piece together the filename (renamed from files provided by Aurina):
    dataFile = sprintf('%s%s%s',dataFileBase,dataFileProbe,dataFileNorm);

    % (1) Load probe information and reformat into a table:
    load(dataFile,'probeInformation');
    entrez_id = probeInformation.EntrezID;
    acronym = probeInformation.GeneSymbol;
    Name = probeInformation.GeneName;
    DSscore = probeInformation.DS;
    geneInfo = table(entrez_id,acronym,Name,DSscore);

    switch gParam.whatParcellation
    case 'APARC'
        load(dataFile,'SampleGeneExpression');
        fprintf(1,'Loaded human expression data from %s\n',which(dataFile));
        whatROI = SampleGeneExpression(:,1);
        sampleExpression = SampleGeneExpression(:,2:end);

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
        if length(ia) < numROIs
            error('ROIs not labeled according to APARC text file IDs...?');
        end
        fprintf(1,'Trimmed structure info to match %u structures\n',...
                                                    height(structInfo));

        fprintf(1,['Matching ROIs to structure info & meaning expression ',...
                        'levels of samples in each ROI\n']);
        numGenes = size(sampleExpression,2);
        geneData = zeros(numROIs,numGenes);
        for i = 1:numROIs
            geneData(i,:) = mean(sampleExpression(whatROI==ROIs(i),:),1);
        end
        fprintf(1,'Done.\n');

    case 'HCP'
        load(dataFile,'geneROI');
        fprintf(1,'Loaded human expression data from %s\n',which(dataFile));
        whatROI = geneROI(:,1);
        geneData = geneROI(:,2:end);
        % Now we need to get the structure info
        ID = whatROI;
        numROIs = length(ID);
        acronym = arrayfun(@(x)sprintf('parcel-%u',ID),ID,'UniformOutput',false);
        % (see Table 1 of the following supplementary info for more information on parcels:
        % https://images-nature-com.ezproxy.lib.monash.edu.au/full/nature-assets/nature/journal/v536/n7615/extref/nature18933-s3.pdf)
        isLeft = true(numROIs,1);
        isCortex = true(numROIs,1);
        structInfo = table(ID,acronym,isLeft,isCortex);
    end
end

%-------------------------------------------------------------------------------
% Further normalization:
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
