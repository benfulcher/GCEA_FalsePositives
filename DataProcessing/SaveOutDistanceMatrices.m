
% Mouse:
distMat = GiveMeDistanceMatrix('mouse','all');
dlmwrite('mouseDistMat.csv',distMat)

% Human cortex:
distMat = GiveMeDistanceMatrix('human');
dlmwrite('humanDistMat.csv',distMat)

%-------------------------------------------------------------------------------
% ----------Run python code from Burt -> mouseSurrogateMaps.csv----------
%-------------------------------------------------------------------------------

VisualizeSurrogateMaps
