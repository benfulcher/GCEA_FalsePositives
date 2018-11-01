function [geneData,geneInfo,structInfo] = LoadMeG(gParam)
% LoadMeG Load gene data as a matrix with tables for metadata info about rows and columns
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs:
if nargin < 1 || isempty(gParam)
    warning('Using default mouse parameters for gene expression')
    humanOrMouse = 'mouse';
    params = GiveMeDefaultParams(humanOrMouse);
    gParam = params.g;
end

%-------------------------------------------------------------------------------
switch gParam.humanOrMouse
case 'surrogate-mouse'
    % dataFileSurrogate = 'mouseSurrogate_rho5_d010.csv';
    dataFileSurrogate = 'mouseSurrogate_rho8_d040.csv';
    dataFileReal = GiveMeFile('AllenMouseGene');

    fprintf(1,'Geometric surrogate mouse maps from %s\n',dataFileSurrogate);
    geneData = dlmread(dataFileSurrogate,',',1,1);
    numFakeGenes = size(geneData,2);

    % Assign random gene metadata:
    fprintf(1,'Real Allen SDK data from %s\n',dataFileReal);
    load(dataFileReal,'geneInfo','structInfo');
    numRealGenes = height(geneInfo);
    rp = sort(randperm(numRealGenes,numFakeGenes));
    fprintf(1,'Assigning a random %u/%u genes\n',numFakeGenes,numRealGenes);
    geneInfo = geneInfo(rp,:);

case 'mouse'
    % Get NEW DATA FROM SDK RETRIEVALS:
    dataFile = GiveMeFile('AllenMouseGene');
    fprintf(1,'New Allen SDK data from %s\n',dataFile);
    load(dataFile,'GeneExpData','geneInfo','structInfo');

    % Use combination (z-score) sections to estimate expression:
    geneData = GeneExpData.combZ.(gParam.energyOrDensity);
case 'human'
    if gParam.normalizeSeparately
        % Cortex and subcortex normalized separately:
        dataFile = '100DS220scaledRobustSigmoidNSGDSQC1LcortexSubcortexSEPARATE_ROI_NOdistCorrSurfaceANDEuclidean.mat';
    else
        % Cortex and subcortex normalized together:
        dataFile = '100DS220scaledRobustSigmoidNSGDSQC1LcortexSubcortex_ROI_NOdistCorrSurfaceANDEuclidean.mat';
    end

    % Get ROI x gene matrix (separate off first colum for the ROI ID)
    load(dataFile,'parcelExpression');
    ROI_ID = parcelExpression(:,1); % ROI IDs (from cust100)
    isCortex = (ROI_ID<=100);
    structInfo = table(ROI_ID,isCortex);
    geneData = parcelExpression(:,2:end);

    % Now we need information about the columns of the gene expression matrix
    % (read in and format as a table):
    load(dataFile,'probeInformation');
    entrez_id = probeInformation.EntrezID;
    acronym = probeInformation.GeneSymbol;
    probeName = probeInformation.ProbeName;
    DS_score = probeInformation.DS;
    geneInfo = table(entrez_id,acronym,probeName,DS_score);

case 'human-old'
    % Data provided by Aurina in 2017
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
% Apply a structure filter
%-------------------------------------------------------------------------------
if ~isempty(structInfo)
    [~,geneData,structInfo] = filterStructures(gParam.structFilter,structInfo,[],geneData);
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
    fprintf(1,'--Normalized expression across each brain region using %s\n',gParam.normalizationRegion);
end

%-------------------------------------------------------------------------------
% Subset genes:
%-------------------------------------------------------------------------------
if ~isempty(gParam.subsetOfGenes)
    warning('Only looking at a random set of %u genes',gParam.subsetOfGenes);
    rp = randperm(size(geneData,2));
    rp = rp(1:gParam.subsetOfGenes);
    geneData = geneData(:,rp);
    geneInfo = geneInfo(rp,:);
end

end
