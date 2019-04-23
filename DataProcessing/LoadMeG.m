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
case 'mouse'
    % Get NEW DATA FROM SDK RETRIEVALS:
    dataFile = GiveMeFile('AllenMouseGene');
    fprintf(1,'New Allen SDK data from %s\n',dataFile);
    load(dataFile,'GeneExpData','geneInfo','structInfo');

    % Use combination (z-score) sections to estimate expression:
    geneData = GeneExpData.combZ.(gParam.energyOrDensity);

case 'human'
    switch gParam.whatParcellation
    case 'HCP'
        dataFile = GiveMeFile('HumanGene_HCP');
    case 'cust100'
        if gParam.normalizeSeparately
            dataFile = GiveMeFile('HumanGene_cust100_normSeparate');
        else
            dataFile = GiveMeFile('HumanGene_cust100_normTogether');
        end
    otherwise
        error('Unknown parcellation: %s',gParam.whatParcellation);
    end

    % Get ROI x gene matrix (separate off first colum for the ROI ID)
    load(dataFile,'parcelExpression');
    ROI_ID = parcelExpression(:,1); % ROI IDs (from cust100)
    switch gParam.whatParcellation
    case 'HCP'
        isCortex = true(size(ROI_ID));
    case 'cust100'
        isCortex = (ROI_ID<=100);
    end
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

case {'surrogate-mouse','surrogate-human'}
    % Start by loading in the real data:
    switch gParam.humanOrMouse
    case 'surrogate-mouse'
        gParam.humanOrMouse = 'mouse';
    case 'surrogate-human'
        gParam.humanOrMouse = 'human';
    end
    [geneDataReal,geneInfo,structInfo] = LoadMeG(gParam);
    numRealGenes = height(geneInfo);
    numAreas = height(structInfo);

    % Surrogate maps:
    switch gParam.whatSurrogate
    case 'spatialLag'
        % Get the pre-computed surrogate data:
        switch gParam.humanOrMouse
        case 'surrogate-mouse'
            dataFileSurrogate = 'mouseSurrogate_N10000_rho8_d040.csv';
        case 'surrogate-human'
            dataFileSurrogate = 'humanSurrogate_N10000_rho8_d02000.csv';
        end
        % Surrogate maps pre-generated using the spatial lag model:
        fprintf(1,'Surrogate brain maps from the spatial lag model, from %s\n',dataFileSurrogate);
        geneData = dlmread(dataFileSurrogate,',',1,1);
        numFakeGenes = size(geneData,2);
        % Assign random gene metadata:
        rp = randperm(numRealGenes,numFakeGenes);
        fprintf(1,'Assigning metadata to genes at RANDOM (%u/%u genes)\n',numFakeGenes,numRealGenes);
        geneInfo = geneInfo(rp,:);

    case 'randomUniform'
        % Uniformly distributed numbers between 0 and 1
        geneData = rand(size(geneDataReal));

    case 'independentSpatialShuffle'
        % Surrogate maps generated through (independent) random shuffling across brain areas
        % (should be consistent with random noise)
        fprintf(1,'Surrogate brain maps from independent random shuffling\n');
        geneData = geneDataReal;
        for j = 1:numRealGenes
            rp = randperm(numAreas);
            geneData(:,j) = geneDataReal(rp,j);
        end

    case 'coordinatedSpatialShuffle'
        % Surrogate maps generated through random shuffling across brain areas
        fprintf(1,'Surrogate brain maps from coordinated random shuffling\n');
        rp = randperm(numAreas);
        geneData = geneDataReal(rp,:);

    case 'geneShuffle'
        % Random shuffling of genes (randomizing association with metadata)
        fprintf(1,'Surrogate brain maps from independent random shuffling\n');
        rp = randperm(numRealGenes);
        geneData = geneDataReal(:,rp);

    otherwise
        error('Unknown surrogate method: ''%s''',whatSurrogate);
    end
end

%-------------------------------------------------------------------------------
% Apply a structure filter
%-------------------------------------------------------------------------------
if ~isempty(structInfo)
    [~,geneData,structInfo] = filterStructures(gParam.structFilter,structInfo,[],geneData);
end

%-------------------------------------------------------------------------------
% Filter on goodness
%-------------------------------------------------------------------------------
minGoodPropGene = 0.3;
minGoodPropArea = 0.1;

% Filter genes on threshold:
nanPropGene = mean(isnan(geneData));
isGoodGene = (nanPropGene < minGoodPropGene);
fprintf(1,'Keeping %u/%u good genes\n',sum(isGoodGene),length(isGoodGene));
geneInfo = geneInfo(isGoodGene,:);
geneData = geneData(:,isGoodGene);

% Filter areas on threshold:
nanPropArea = mean(isnan(geneData),2);
isGoodArea = (nanPropArea < minGoodPropArea);
fprintf(1,'Keeping %u/%u good areas\n',sum(isGoodArea),length(isGoodArea));
structInfo = structInfo(isGoodArea,:);
geneData = geneData(isGoodArea,:);

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
