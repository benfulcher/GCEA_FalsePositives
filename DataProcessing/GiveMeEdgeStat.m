function [edgeData,regionInfo] = GiveMeEdgeStat(connectomeType,pThreshold,whatEdgeProperty,correctDistance,numThresholds)
% Returns edge data for a given connectome and edge property
%-------------------------------------------------------------------------------

if nargin < 4
    correctDistance = false;
end
if nargin < 5
    numThresholds = 10;
end

%-------------------------------------------------------------------------------
% First get the adjacency matrix for the connectome shorthand specified:
switch connectomeType
case 'Oh-brain'
    [A_bin,regionInfo] = GiveMeAdj('Oh',pThreshold,true,'right',false);
    A_wei = GiveMeAdj('Oh',pThreshold,false,'right',false);
case 'Oh-cortex'
    [A_bin,regionInfo] = GiveMeAdj('Oh',pThreshold,true,'right',true);
    A_wei = GiveMeAdj('Oh',pThreshold,false,'right',true);
otherwise
    error('Unknown connectome: %s',connectomeTypes);
end

%-------------------------------------------------------------------------------
% Compute an edge-level statistic from this connectome:
switch whatEdgeProperty
case 'wei-communicability'
    edgeData = communicability(A_wei);
    edgeData(A_wei==0) = 0; % only put on real edges
case 'bin-communicability'
    edgeData = communicability(A_bin);
    edgeData(~A_bin) = 0; % only put on real edges
case 'bin-betweenness'
    edgeData = edge_betweenness_bin(A_bin);
case 'wei-betweenness'
    edgeData = edge_betweenness_wei(A_bin);
case 'distance'
    C = load('Mouse_Connectivity_Data.mat','Dist_Matrix');
    edgeData = C.Dist_Matrix{1,1}/1000; % ipsilateral distances in the right hemisphere
    edgeData(tril(true(size(dData)),-1)) = 0; % remove lower diagonal (symmetric)
otherwise
    error('Unknown edge property: %s',whatEdgeProperty);
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
