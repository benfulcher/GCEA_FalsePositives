function [regionStruct,connProp] = GiveMeConnectivityData(whatConnectivityData,regionFilter)
% Idea is to see which connectivity properties are related to cell type
% differences
%-------------------------------------------------------------------------------
% Connectivity data options:
% 'Ypma-complete-ipsi': subset of complete ipsilateral cortical connectome
%                       estimated by Ypma & Bullmore
% 'Oh-ipsi': ipsilateral connectivity from Oh et al.
%
% Will restrict the connectome to the regions given in regionSubset if provided

%-------------------------------------------------------------------------------
% Parameters
if nargin < 1
    whatConnectivityData = 'Oh-ipsi'; % isocortex from Oh et al.
end
if nargin < 2
    regionFilter = [];
end

%-------------------------------------------------------------------------------
% First get the data:
switch whatConnectivityData
case 'Ypma-ipsi'
    [W_rect,sourceRegions,targetRegions] = ImportCorticalConnectivityWeights();
    [W,regionNames] = MakeCompleteConnectome(W_rect,sourceRegions,targetRegions);
    [regionStruct,ia] = MatchRegionsOh([],regionNames);
    W = W(ia,ia);
case 'Oh-ipsi'
    C = load('Mouse_Connectivity_Data.mat');
    W = GiveMeAdj(C,'zero','ipsi',0,0.05);
    regionStruct = C.RegionStruct;
end

%-------------------------------------------------------------------------------
% Match to a regional subset:
%-------------------------------------------------------------------------------
[regionStruct,ix] = StructureFilter(regionStruct,regionFilter);
W = W(ix,ix);

%-------------------------------------------------------------------------------
% Make a zeroed version (NaNs -> 0):
Wzero = W;
Wzero(isnan(Wzero)) = 0;
% Make a binary zeroed version:
% W_binaryZero = Wzero;
% W_binaryZero(W_binaryZero > 0) = 1;

% Compute conectivity properties:
connProp = table();
connProp.inStrength = sum(Wzero,1)';
connProp.outStrength = sum(Wzero,2);
connProp.strength = sum(Wzero,1)' + sum(Wzero,2);
% connProp.inDegree = sum(W_binaryZero,1)';
% connProp.outDegree = sum(W_binaryZero,2);

end
