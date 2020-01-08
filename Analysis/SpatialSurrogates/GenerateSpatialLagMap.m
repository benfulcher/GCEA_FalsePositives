function X = GenerateSpatialLagMap(distMat,d0,rho,numMaps)

% Set defaults
if nargin < 3 || isempty(rho)
    rho = 0.8;
end
if nargin < 4
    numMaps = 1;
end
%-------------------------------------------------------------------------------

numPoints = size(distMat,1);

% W_{ij} = exp(-d_{ij} / d0)
W = exp(-distMat / d0);

% Clear diagonal:
identity = logical(eye(numPoints));
W(identity) = 0;

% Random numbers:
X = zeros(numPoints,numMaps);

if numMaps > 1
    fprintf(1,'Generating %u (null) spatially autocorrelated maps\n',numMaps);
end

for i = 1:numMaps
    u = randn(numPoints,1);
    % x = (identity - rho * W)\u;
    % x = inv(identity - rho * W)*u;
    X(:,i) = rho * W * u;
end

end
