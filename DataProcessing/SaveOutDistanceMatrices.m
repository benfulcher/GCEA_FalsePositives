
% Mouse:
distMat = GiveMeDistanceMatrix('mouse','all');
dlmwrite('mouseDistMat.csv',distMat)

% Mouse:
distMat = GiveMeDistanceMatrix('mouse','cortex');
dlmwrite('mouseCortexDistMat.csv',distMat)

% Human cortex:
distMat = GiveMeDistanceMatrix('human');
numAreas = length(distMat);
dlmwrite(sprintf('humanDistMat_%u.csv',numAreas),distMat)

%-------------------------------------------------------------------------------
% ----------Run python code from Burt -> mouseSurrogateMaps.csv----------
%-------------------------------------------------------------------------------

% VisualizeSurrogateMaps
