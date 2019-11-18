%-------------------------------------------------------------------------------
% Pairwise Euclidean distance matrices
%-------------------------------------------------------------------------------

% Mouse brain:
distMat = GiveMeDistanceMatrix('mouse','all');
dlmwrite('mouseDistMat.csv',distMat)

% Mouse cortex:
distMat = GiveMeDistanceMatrix('mouse','cortex');
dlmwrite('mouseCortexDistMat.csv',distMat)

% Human cortex:
distMat = GiveMeDistanceMatrix('human');
numAreas = length(distMat);
dlmwrite(sprintf('humanDistMat_%u.csv',numAreas),distMat)

%-------------------------------------------------------------------------------
% Degree phenotype:
%-------------------------------------------------------------------------------
doBinarize = true;

% Mouse brain:
params = GiveMeDefaultParams('mouse','all');
[k,structInfoConn] = ComputeDegree(params,doBinarize);
dlmwrite('mouseDegree.csv',k)

% Mouse cortex:
params = GiveMeDefaultParams('mouse','cortex');
[k,structInfoConn] = ComputeDegree(params,doBinarize);
dlmwrite('mouseCortexDegree.csv',k)

% Human cortex:
params = GiveMeDefaultParams('human','cortex');
[k,structInfoConn] = ComputeDegree(params,doBinarize);
numAreas = length(k);
dlmwrite(sprintf('humanDegree_%u.csv',numAreas),k)

%-------------------------------------------------------------------------------
% Some internal analyses
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% ----------Run python code from Burt -> mouseSurrogateMaps.csv----------
%-------------------------------------------------------------------------------

% VisualizeSurrogateMaps
