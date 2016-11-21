function theAdjMat = GiveMeAdj(C,whatAdj,whatHemisphere,justCortex,pThreshold)
% Gives a string identifying the type of normalization to apply, then returns
% the gene data for that normalization.
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-07-17
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Check Inputs
% ------------------------------------------------------------------------------

if nargin < 3 || isempty(whatHemisphere)
    whatHemisphere = 'ipsi';
end

if nargin < 4 || isempty(justCortex)
    justCortex = 0;
end

if nargin < 5 || isempty(pThreshold)
    pThreshold = 0.05;
end

% ------------------------------------------------------------------------------
% Take the appropriate hemisphere
% ------------------------------------------------------------------------------

switch whatHemisphere
case 'ipsi'
    theIndex = 1;
case 'contra'
    theIndex = 2;
end

% ------------------------------------------------------------------------------
% Retrieve and process the data
% ------------------------------------------------------------------------------
switch whatAdj
case 'raw'
    % No normalization, or p-filtering
    theAdjMat = C.Conn_W{theIndex};

case 'zero'
    % Raw weights but with zeros where either zero or NaN
    theAdjMat = C.Conn_W{theIndex};

    % Remove diagonal entries:
    theAdjMat(logical(eye(size(theAdjMat)))) = 0;

    % Zero high-p links using the given p-threshold:
    theAdjMat = filterP(theAdjMat);

case 'binary'
    % Raw weights but with zeros where either zero or NaN
    theAdjMat = C.Conn_W{theIndex};

    % Remove diagonal entries:
    theAdjMat(logical(eye(size(theAdjMat)))) = 0;

    % Zeros for 0.05-significant connections, no diagonal
    theAdjMat = filterP(theAdjMat);

    % Make binary:
    theAdjMat(theAdjMat>0) = 1;

otherwise
    error('Unknown adjacency matrix option: ''%s''',whatAdj);
end

% ------------------------------------------------------------------------------
% Take just cortex?
% ------------------------------------------------------------------------------
if justCortex
    isCortex = strcmp({C.RegionStruct.MajorRegionName},'Isocortex');
    theAdjMat = theAdjMat(isCortex,isCortex);
end

% ------------------------------------------------------------------------------
function AdjOut = filterP(AdjIn)
    % Sets high-p-value links to zero:
    AdjOut = AdjIn;
    AdjOut(C.Conn_p{theIndex} > pThreshold) = 0;
end

end
