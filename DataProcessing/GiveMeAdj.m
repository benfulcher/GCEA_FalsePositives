function [theAdjMat,regionAcronyms,adjPVals] = GiveMeAdj(whatData,pThreshold,doBinarize,...
                                    whatWeightMeasure,whatHemispheres,whatFilter)
% Gives a string identifying the type of normalization to apply, then returns
% the gene data for that normalization.
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check Inputs
%-------------------------------------------------------------------------------
if nargin < 1
    whatData = 'mouse-Oh';
    % whatData = 'human-HCP';
end
if nargin < 2 || isempty(pThreshold)
    pThreshold = 0.05;
end
if nargin < 3
    doBinarize = false;
end
if nargin < 4
    whatWeightMeasure = 'NCD'; % normalized connection density
end
if nargin < 5
    whatHemispheres = 'right';
end
if nargin < 6
    whatFilter = 'all';
end

%-------------------------------------------------------------------------------
% Load in and minimally preprocess the data:
if isstruct(whatData)
    C = whatData;
    whatData = 'mouse-Oh';
end
switch whatData
case 'mouse-Oh'
    if ~exist('C','var')
        C = load('Mouse_Connectivity_Data.mat','Conn_W','Conn_p','regionAcronyms');
    end
    % Ipsi:
    switch whatHemispheres
    case 'right'
        ind = [1,1];
    case 'left'
        ind = [2,2];
    end
    theAdjMat = C.Conn_W{ind(1),ind(2)};
    adjPVals = C.Conn_p{ind(1),ind(2)};
    % Remove diagonal entries:
    theAdjMat(logical(eye(size(theAdjMat)))) = 0;
    % Zero high-p links using the given p-threshold:
    theAdjMat = filterP(theAdjMat,adjPVals);

    % Set custom edge weight measure:
    numROIs = length(theAdjMat);
    switch whatWeightMeasure
    case 'NCS' % normalised connection strength
        fprintf(1,'~~Normalized connection strength~~\n');
        theAdjMat = theAdjMat;
    case 'NCD' % normalized connection density
        fprintf(1,'~~Normalized connection density~~\n');
        % Divide by destination ROI size, Y
        roi_volume = GetROIVolumes(C);
        for i_target = 1:numROIs
            theAdjMat(:,i_target) = theAdjMat(:,i_target)/roi_volume(i_target);
        end
    case 'CS' % connection strength
        fprintf(1,'~~Connection strength~~\n');
        % Multiply by source (row) by ROI volume
        roi_volume = GetROIVolumes(C);
        numROIs = length(theAdjMat);
        for i_source = 1:numROIs
            theAdjMat(i_source,:) = theAdjMat(i_source,:)*roi_volume(i_source);
        end
    case 'CD' % connection density
        fprintf(1,'~~Connection density~~\n');
        % multiply each weight by the source volume and divide by target volume
        roi_volume = GetROIVolumes(C);
        numROIs = length(theAdjMat);
        % multiply by source volume:
        for i_source = 1:numROIs
            theAdjMat(i_source,:) = theAdjMat(i_source,:)*roi_volume(i_source);
        end
        % divide by target volume:
        for i_target = 1:numROIs
            theAdjMat(:,i_target) = theAdjMat(:,i_target)/roi_volume(i_target);
        end
    otherwise
        error('Unknown edge weight type: ''%s''',whatWeights);
    end

    regionAcronyms = C.regionAcronyms;

case 'mouse-Ypma'
    [W_rect,sourceRegions,targetRegions] = ImportCorticalConnectivityWeights();
    [W,regionNames] = MakeCompleteConnectome(W_rect,sourceRegions,targetRegions);
    [regionStruct,ia] = MatchRegionsOh([],regionNames);
    W = W(ia,ia);
    theAdjMat = W;
    adjPVals = [];

case {'human-HCP-HCP','human-HCP-APARC'}
    % Get the 100-HCP group connectome [Group connectomes are made by removing least
    % consistent (in terms of weight) links to reach the average density of
    % individual connectomes.]
    switch whatData
    case 'human-HCP-HCP'
        % (HCP parcellation)
        fprintf(1,'HCP parcellation\n');
        C = load('HCPgroupConnectomesHCP.mat');
        structInfo = GiveMeHCPNames();
    case 'human-HCP-APARC'
        % (APARC parcellation)
        fprintf(1,'APARC parcellation\n');
        C = load('HCPgroupConnectomesAparcaseg.mat');
        structInfo = GiveMeAPARCNames();
    end

    % Extract the right weight measure:
    switch  whatWeightMeasure
    case 'count'
        theAdjMat = C.AdjCount;
    case 'density'
        theAdjMat = C.AdjDens;
    otherwise
        error('Unknown weight measure: ''%s''',whatWeightMeasure);
    end

    % Take just left hemisphere:
    if strcmp(whatData,'human-HCP-HCP')
        % Left cortex is 1:180 (180)
        % Left subcortex 181:190 (10)
        % Right cortex 191:370 (180)
        % Right subcortex 371:380 (10)
        isLeftCortex = 1:180;
        theAdjMat = theAdjMat(isLeftCortex,isLeftCortex);
    end

    adjPVals = []; % Just to fill this output (not used for human data)

    % Get ROI names (regionAcronyms):
    regionAcronyms = structInfo.acronym;

otherwise
    error('Unknown data source, ''%s''',whatData);
end

%-------------------------------------------------------------------------------
% Binarize
if doBinarize
    theAdjMat = theAdjMat;
    theAdjMat(theAdjMat > 0) = 1;
end

%-------------------------------------------------------------------------------
% Filter structures
if ismember(whatFilter,{'isocortex','cortex'})
    if strcmp(whatData(1:5),'mouse')
        % Need to load in structInfo (pretty inefficient -- throw out large loaded gene data)
        [~,~,structInfo] = LoadMeG();
        % Match to regions:
        keepStruct = strcmp(structInfo.divisionLabel,'Isocortex');
        theAdjMat = theAdjMat(keepStruct,keepStruct);
    elseif strcmp(whatData(1:5),'human')
        isCTX = regionStruct.isCortex;
        theAdjMat = theAdjMat(isCTX,isCTX);
        regionStruct = regionStruct(isCTX,:);
        regionAcronyms = regionAcronyms(isCTX);
    else
        error('Unknown label for organism -- are you a human or are you a mouse?');
    end
end

if strcmp(whatData(1:5),'human')
    switch whatHemispheres
    case 'right'
        % Keep right hemisphere only (human data)
        isRight = structInfo.isRight;
        theAdjMat = theAdjMat(isRight,isRight);
        regionAcronyms = regionAcronyms(isRight);
        structInfo = structInfo(isRight,:);
        fprintf(1,'Filtered to %u right-hemisphere ROIs\n',sum(isRight));
    case 'left'
        % Keep left hemisphere only (human data)
        isLeft = structInfo.isLeft;
        theAdjMat = theAdjMat(isLeft,isLeft);
        regionAcronyms = regionAcronyms(isLeft);
        structInfo = structInfo(isLeft,:);
        fprintf(1,'Filtered to %u left-hemisphere ROIs\n',sum(isLeft));
    case 'both'
        % yeah
    otherwise
        error('What do you mean %s hemisphere?!',whatHemispheres);
    end
end

% ------------------------------------------------------------------------------
function AdjThresh = filterP(AdjIn,pValues)
    % Sets high-p-value links to zero:
    AdjThresh = AdjIn;
    AdjThresh(pValues > pThreshold) = 0;
end

end
