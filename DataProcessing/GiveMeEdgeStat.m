function [edgeData,regionAcronyms] = GiveMeEdgeStat(connectomeType,pThreshold,whatWeightMeasure,...
                                whatEdgeProperty,correctDistance,numThresholds)
% Returns edge data for a given connectome and edge property
% [[REDUNDANT FOR GiveMeEdgeMeasure]]
%-------------------------------------------------------------------------------
if nargin < 2
    pThreshold = 0.05;
end
if nargin < 3
    whatWeightMeasure = 'NCD';
end
if nargin < 5
    correctDistance = false;
end
if nargin < 6
    numThresholds = 10;
end

%-------------------------------------------------------------------------------
% First get the adjacency matrix for the connectome shorthand specified:
switch connectomeType
case 'Oh-brain'
    [A_wei,regionAcronyms] = GiveMeAdj('Oh',pThreshold,false,whatWeightMeasure,'right','all');
    A_bin = double(A_wei > 0);
case 'Oh-cortex'
    [A_wei,regionAcronyms] = GiveMeAdj('Oh',pThreshold,false,whatWeightMeasure,'right','cortex');
    A_bin = double(A_wei > 0);
otherwise
    error('Unknown connectome: %s',connectomeTypes);
end

%-------------------------------------------------------------------------------
% Compute an edge-level statistic from this connectome:
onlyOnEdges = true;
A_p = [];
if strcmp(whatEdgeProperty,'distance')
    C = load('Mouse_Connectivity_Data.mat','Dist_Matrix');
    edgeData = C.Dist_Matrix{1,1}/1000; % ipsilateral distances in the right hemisphere
    edgeData(tril(true(size(dData)))) = 0; % remove lower diagonal (symmetric)
else
    edgeData = GiveMeEdgeMeasure(whatEdgeProperty,A_bin,A_wei,onlyOnEdges);
end

%-------------------------------------------------------------------------------
% Distance correction:
% (based on method in BF_PlotQuantiles)
if correctDistance
    connValues = (edgeData > 0); % connections exist here
    C = load('Mouse_Connectivity_Data.mat','Dist_Matrix');
    xData = C.Dist_Matrix{1,1}(connValues)/1000; % distances for edges that exist
    yData = edgeData(connValues); % edge statistic for edges that exist

    xThresholds = arrayfun(@(x)quantile(xData,x),linspace(0,1,numThresholds));
    xThresholds(end) = xThresholds(end)*1.1; % make sure all data included in final bin

    yDataCorrected = nan(size(yData));
    for p = 1:numThresholds-1
        inBin = (xData>=xThresholds(p) & xData < xThresholds(p+1));
        yDataCorrected(inBin) = yData(inBin) - mean(yData(inBin));
    end
    if any(isnan(yDataCorrected))
        keyboard
        error('NaNs in edge data :/');
    end

    edgeDataCorrected = zeros(size(edgeData));
    edgeDataCorrected(connValues) = yDataCorrected;
end

end
