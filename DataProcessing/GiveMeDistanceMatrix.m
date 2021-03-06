function distanceMat = GiveMeDistanceMatrix(params)
% Gives a pairwise distance matrix
%-------------------------------------------------------------------------------
if nargin < 1
    params = GiveMeDefaultParams();
end

%-------------------------------------------------------------------------------
switch params.humanOrMouse
case 'mouse'
    % Pairwise ipsilateral Euclidean distance data:
    % warning('Distances only work at the moment for whole-brain data')
    load('Mouse_Connectivity_Data.mat','Dist_Matrix','regionAcronyms');
    distanceMat = Dist_Matrix{1,1}/1000;
    acronym = regionAcronyms;

    % Map to division labels:
    dataFile = GiveMeFile('AllenMouseGene');
    load(dataFile,'structInfo');

case 'human'
    % Latest data processed by Aurina
    % cust100: cortical and subcortical structures normalized separately:
    switch params.g.whatParcellation
    case 'HCP'
        dataFile = GiveMeFile('HumanGene_HCP');
    case 'cust100'
        if params.g.normalizeSeparately
            dataFile = GiveMeFile('HumanGene_cust100_normSeparate');
        else
            dataFile = GiveMeFile('HumanGene_cust100_normTogether');
        end
    end
    load(dataFile,'averageDistance');
    distanceMat = (averageDistance + averageDistance')/2;
    ROIs_distance = 1:length(distanceMat);

    % Get gene-expression data to ensure a match:
    [geneData,geneInfo,structInfo] = LoadMeG(params.g);

    % Get gene-expression data to ensure a match:
    fprintf(1,'Keeping distances to correspond to the brain areas kept for the default gene-expression data\n');
    doKeep = ismember(ROIs_distance,structInfo.ROI_ID);
    warning('Keeping %u/%u regions for distance information',sum(doKeep),length(doKeep));
    distanceMat = distanceMat(doKeep,doKeep);
    params.structFilter = 'all';

    % Keep only cortex:
    % warning('Hard keeping first 100 cust100 cortical parcels (out of %u)',length(distanceMat));
    % distanceMat = distanceMat(1:100,1:100);

case 'human2017'
    % Results using data provided by Aurina in 2017
    warning('ASSUMING HCP LEFT CORTICAL PARCELLATION!!! :-O')
    fileName = '360parcellationLcortex_ProbeMean.mat';
    fprintf(1,'Loading data from ''%s''\n',fileName);
    load(fileName,'coordinatesROI');
    fprintf(1,'%u ROIs\n',size(coordinatesROI,1));
    ROIIDs = coordinatesROI(:,1);
    % Check that they match the expression data
    load(fileName,'geneROI');
    expression_ROI_IDs = geneROI(:,1);
    if ~all(ROIIDs==expression_ROI_IDs)
        error('Error matching ROI IDs for %s',fileName);
    end
    coOrds = coordinatesROI(:,2:end);
    distanceMat = squareform(pdist(coOrds,'Euclidean'));
end

%-------------------------------------------------------------------------------
% Filter structures:
if ~strcmp(params.structFilter,'all')
    fprintf(1,'Filtering structures...?\n');
    [~,~,~,keepStruct] = filterStructures(params.structFilter,structInfo);
    distanceMat = distanceMat(keepStruct,keepStruct);
end

%-------------------------------------------------------------------------------

fprintf(1,'%u x %u pairwise distance matrix\n',size(distanceMat,1),size(distanceMat,2));

end
