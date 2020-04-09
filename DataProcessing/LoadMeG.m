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
        paramsFull = GiveMeDefaultParams(gParam.humanOrMouse,gParam.structFilter);
        dataFileSurrogate = FindNullFile(paramsFull);
        % Surrogate maps pre-generated using the spatial lag model:
        fprintf(1,'SUBSTITUTING WITH Surrogate brain maps from the spatial lag model, from %s\n',...
                            dataFileSurrogate);
        load(dataFileSurrogate,'nullMaps');
        geneData = nullMaps;
        % Reshape to match any filtering applied to the loaded data:
        if size(geneData,1) > size(geneDataReal,1)
            error('Data do not match :-/')
        end

        % Assign random gene metadata:
        numFakeGenes = size(geneData,2);
        if numFakeGenes > numRealGenes
            geneData = geneData(:,1:numRealGenes);
        else
            rp = randperm(numRealGenes,numFakeGenes);
            fprintf(1,'Assigning metadata to genes at RANDOM (%u/%u genes)\n',numFakeGenes,numRealGenes);
            geneInfo = geneInfo(rp,:);
        end

    case 'randomMap'
        geneData = ShuffleMyMatrix(geneDataReal,'randomUniform');

    case 'independentSpatialShuffle'
        % Surrogate maps generated through (independent) random shuffling across brain areas
        % (should be consistent with random noise)
        fprintf(1,'Surrogate brain maps from independent random shuffling\n');
        geneData = ShuffleMyMatrix(geneDataReal,'independentRowShuffle');

    case 'coordinatedSpatialShuffle'
        % Surrogate maps generated through random shuffling across brain areas
        fprintf(1,'Surrogate brain maps from coordinated random shuffling\n');
        geneData = ShuffleMyMatrix(geneDataReal,'coordinatedRowShuffle');

    case 'geneShuffle'
        % Random shuffling of genes (randomizing association with metadata)
        fprintf(1,'Surrogate brain maps from independent random shuffling\n');
        geneData = ShuffleMyMatrix(geneDataReal,'coordinatedColumnShuffle');

    otherwise
        error('Unknown surrogate method: ''%s''',gParam.whatSurrogate);
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
% Filter genes on threshold:
notNanPropGene = mean(~isnan(geneData));
isGoodGene = (notNanPropGene >= gParam.minGoodPropGene);
fprintf(1,'Keeping %u/%u good genes\n',sum(isGoodGene),length(isGoodGene));
geneInfo = geneInfo(isGoodGene,:);
geneData = geneData(:,isGoodGene);

% Filter areas on threshold:
notNanPropArea = mean(~isnan(geneData),2);
isGoodArea = (notNanPropArea >= gParam.minGoodPropArea);
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
