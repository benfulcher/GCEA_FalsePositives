function distanceMat = GiveMeDistanceMatrix(whatSpecies)
% Gives a pairwise distance matrix
%-------------------------------------------------------------------------------

switch whatSpecies
case 'mouse'
    % Pairwise ipsilateral Euclidean distance data:
    % warning('Distances only work at the moment for whole-brain data')
    % (could add filtering)
    load('Mouse_Connectivity_Data.mat','Dist_Matrix');
    distanceMat = Dist_Matrix{1,1}/1000;
case 'human'
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

fprintf(1,'%u x %u pairwise distance matrix\n',...
                size(distanceMat,1),size(distanceMat,2));

end
