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
distMat = GiveMeDistanceMatrix('mouse','all');
params = GiveMeDefaultParams('mouse','all');
[k,structInfoConn] = ComputeDegree(params,doBinarize);
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
k_bc = boxcox(k);
k_z = zscore(k_bc);
K_Z = k_z*k_z';
% Difference (rather than product-against-mean) approach:
% K_Z_diff = zeros(length(k));
% for i = 1:length(k)
%     for j = 1:length(k)
%         K_Z_diff(i,j) = abs(k(i)-k(j));
%     end
% end
% upperMask = triu(true(size(K_Z)),+1);

f = figure('color','w');
plot(distMat(upperMask),K_Z(upperMask),'.k')
xlabel('Distance')
ylabel('Degree Similarity')

%-------------------------------------------------------------------------------
options = optimset('Display','iter');
SpatialLagHere = @(rho_d0) SpatialLagForm(k_z,distMat,rho_d0(1),rho_d0(2),true);
rho_d0_init = [0.4,0.5];
[opt_param,d] = fminsearch(SpatialLagHere,rho_d0_init,options);
rho = opt_param(1);
d0 = opt_param(2);

options = optimset('Display','iter');
x0 = [0.4,0.5];

k = rho*

%-------------------------------------------------------------------------------
% ----------Run python code from Burt -> mouseSurrogateMaps.csv----------
%-------------------------------------------------------------------------------

% VisualizeSurrogateMaps
