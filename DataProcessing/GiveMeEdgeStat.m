function edgeData = GiveMeEdgeStat(connectomeType,pThreshold,whatEdgeProperty,correctDistance,numThresholds)
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
    A_bin = GiveMeAdj('Oh',pThreshold,true,'right',false);
    A_wei = GiveMeAdj('Oh',pThreshold,false,'right',false);
case 'Oh-cortex'
    A_bin = GiveMeAdj('Oh',pThreshold,true,'right',true);
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
    edgeData = C.Dist_Matrix{1,1}/1000; % ipsilateral distances in the right hemisphere
end

%-------------------------------------------------------------------------------
% Distance correction:
% (based on method in BF_PlotQuantiles)
if correctDistance
    connValues = (edgeData > 0); % connections exist here
    xData = C.Dist_Matrix{1,1}(connValues)/1000;
    yData = edgeData(connValues);

    xThresholds = arrayfun(@(x)quantile(xData,x),linspace(0,1,numThresholds));
    xThresholds(end) = xThresholds(end) + eps; % make sure all data included in final bin

    edgeDataCorrected = nan(size(yData));
    for p = 1:numThresholds-1
        inBin = (xData>=xThresholds(p) & xData < xThresholds(p+1));
        edgeDataCorrected(inBin) = yData(inBin) - mean(yData(inBin));
    end

    edgeData = edgeDataCorrected;
end

end
