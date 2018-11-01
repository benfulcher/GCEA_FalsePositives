
distMat = GiveMeDistanceMatrix('mouse','all');
dlmwrite('mouseDistMat.csv',distMat)

%-------------------------------------------------------------------------------
% ----------Run python code from Burt -> mouseSurrogateMaps.csv----------
%-------------------------------------------------------------------------------

VisualizeSurrogateMaps
