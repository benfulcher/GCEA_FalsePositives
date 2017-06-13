function [edgeData,edgeDataCorrected] = DistanceDependenceEdge(whatEdgeProperty,pThreshold,connectomeType,plotFull)
%-------------------------------------------------------------------------------
% Aim is to plot the distance dependence of edge metrics
%-------------------------------------------------------------------------------

if nargin < 1
    whatEdgeProperty = 'bin_communicability';
end
if nargin < 2
    pThreshold = 0.05;
end
if nargin < 3
    connectomeType = 'Oh-brain';
end
if nargin < 4
    plotFull = true;
end
numQuantiles = 11;

%===============================================================================
C = load('Mouse_Connectivity_Data.mat','Dist_Matrix');
d = C.Dist_Matrix{1,1}/1000; % ipsilateral distances in the right hemisphere

A_bin = GiveMeAdj('Oh',pThreshold,true,'NCD','right');
A_wei = GiveMeAdj('Oh',pThreshold,false,'NCD','right');

switch connectomeType
case 'Oh-brain'
    % nothing to do
case 'Oh-cortex'
    [~,~,structInfo] = LoadMeG({'none','none'},'energy');
    keepStruct = strcmp(structInfo.divisionLabel,'Isocortex');
    structInfo = structInfo(keepStruct,:);
    A_bin = A_bin(keepStruct,keepStruct);
    d = d(keepStruct,keepStruct);
otherwise
    error('Unknown connectome type: %s',connectomeType);
end

% Compute the edge data:
onlyOnEdges = true;
if strcmp(whatEdgeProperty,'distance')
    edgeData = d; % ipsilateral distances in the right hemisphere
else
    edgeData = GiveMeEdgeMeasure(whatEdgeProperty,A_bin,A_wei,onlyOnEdges);
end
connValues = (edgeData ~= 0);

%---------------------------------------------------------------------------
f = figure('color','w');
if plotFull
    subplot(2,2,1)
    histogram(edgeData(connValues));
    xlabel(whatEdgeProperty,'interpreter','none')
    subplot(2,2,2)
    plot(d(connValues),edgeData(connValues),'.k')
    xlabel('d')
    ylabel(whatEdgeProperty,'interpreter','none')
    title(pThreshold)
    subplot(2,2,3:4)
end
edgeDataCorrected = BF_PlotQuantiles(d(connValues),edgeData(connValues),numQuantiles,false,false);
title(sprintf('%s on %u edges',whatEdgeProperty,sum(connValues(:))),'interpreter','none')
xlabel('d');
ylabel(whatEdgeProperty,'interpreter','none')

end
