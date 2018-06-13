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
    warning('Using default mouse parameters for gene expression')
    humanOrMouse = 'mouse';
    params = GiveMeDefaultParams(humanOrMouse);
    gParam = params.g;
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Get NEW DATA FROM SDK RETRIEVALS:
%-------------------------------------------------------------------------------
% try
switch gParam.humanOrMouse
case 'mouse'
    dataFile = '/Users/benfulcher/DropboxSydneyUni/CurrentProjects/CellTypesMouse/Code/Data/AllenGeneDataset_19419.mat';
    fprintf(1,'New Allen SDK-data from %s\n',dataFile);
    load(dataFile,'GeneExpData','geneInfo','structInfo');

    % Use combination (z-score) sections to estimate expression:
    geneData = GeneExpData.combZ.(gParam.energyOrDensity);
case 'human'
    % Piece together filename from parameters:
    switch gParam.whatParcellation
    case 'APARC'
        dataFileBase = 'aparcaseg'; % APARC parcellation
    case 'HCP'
        dataFileBase = '360parcellationLcortex'; % HCP parcellation
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
        fprintf(1,'Loaded sample-wise APARC human expression data from %s\n',which(dataFile));
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
        fprintf(1,'Loaded ROI-wise HCP human expression data from %s\n',which(dataFile));
        whatROI = geneROI(:,1);
        geneData = geneROI(:,2:end);
        % Now we need to get the structure info
        structInfo = GiveMeHCPNames();
        % And filter according to the expression data that is here
        [~,ia,ib] = intersect(whatROI,structInfo.ID,'stable');
        if ~all(diff(ia)==1)
            error('Error matching HCP ROIs');
        end
        structInfo = structInfo(ib,:);
    end
end

%-------------------------------------------------------------------------------
% Further normalization:
%-------------------------------------------------------------------------------
if ~strcmp(gParam.normalizationGene,'none')
    geneData = BF_NormalizeMatrix(geneData,gParam.normalizationGene);
    fprintf(1,'--Normalized expression for each gene using %s\n',gParam.normalizationGene);
end
if ~strcmp(gParam.normalizationRegion,'none')
    geneData = BF_NormalizeMatrix(geneData',gParam.normalizationRegion)';
    fprintf(1,'2. Normalized expression across each brain region using %s\n',gParam.normalizationRegion);
end

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
