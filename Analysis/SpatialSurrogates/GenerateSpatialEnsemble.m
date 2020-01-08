function GenerateSpatialEnsemble(mouseOrHuman,structFilter,doPlot,numMaps)
%===============================================================================
% Matlab Code to generate an ensemble of spatially autocorrelated maps
% (cf. python pipeline in SaveOutDistanceMatrices)
%===============================================================================

%-------------------------------------------------------------------------------
% Check inputs and set defaults:
if nargin < 1
    mouseOrHuman = 'mouse';
end
if nargin < 2
    structFilter = 'all';
end
if nargin < 3
    doPlot = true;
end
if nargin < 4
    numMaps = 1000;
end
%-------------------------------------------------------------------------------

% Get all pairwise distances:
distMat = GiveMeDistanceMatrix(mouseOrHuman,structFilter);
numAreas = length(distMat);

% Check the spatial scale:
params = GiveMeDefaultParams(mouseOrHuman,structFilter);
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
% This normalization by region and area makes the exponential dependence clearer:
geneDataZ = BF_NormalizeMatrix(geneData,'zscore');
geneDataZZ = BF_NormalizeMatrix(geneDataZ','zscore')';
G = corr(geneDataZZ','rows','pairwise');
upperMask = triu(true(size(G)),+1);
xData = distMat(upperMask);
yData = G(upperMask);

% Fit:
[f_handle,Stats,c] = GiveMeFit(xData,yData,'exp0',false);
rho = c.A; % strength of relationship
d0 = 1/c.n; % spatial scale of transcriptional autocorrelation
if doPlot
    % Plot the data and fit:
    f = figure('color','w'); hold('on');
    plot(distMat(upperMask),G(upperMask),'.k')
    xRange = linspace(min(distMat(:)),max(distMat(:)),50);
    plot(xRange,f_handle(xRange),'b','LineWidth',3)
end

fprintf(1,'Mouse-brain gene expression has strength, rho = %g \n',rho);
fprintf(1,'Mouse-brain gene expression has spatial scale, d0 = %g \n',d0);

% Generate null maps with these parameters:
nullMaps = GenerateSpatialLagMap(distMat,d0,rho,numMaps);
if doPlot
    % Compute 2-d projection of data:
    coOrdsXY = mdscale(distMat,2);
    % Plot each surrogate spatial map:
    f = figure('color','w');
    for i = 1:6
        subplot(2,3,i)
        mapNorm = zscore(nullMaps(:,i));
        scatter(coOrdsXY(:,1),coOrdsXY(:,2),25,mapNorm,'filled')
        colormap([flipud(BF_getcmap('blues',9));1,1,1;BF_getcmap('reds',9)])
        xlabel('spatialAxis1')
        ylabel('spatialAxis2')
        title(sprintf('rho = %.3f, d0 = %.3f',rho,d0))
    end
end

%-------------------------------------------------------------------------------
% Save to file
if strcmp(mouseOrHuman,'mouse')
    fileOut = sprintf('%s_%s_Surrogate_N%u.mat',mouseOrHuman,structFilter,numMaps);
else
    fileOut = sprintf('%s_Surrogate_N%u.mat',mouseOrHuman,structFilter,numMaps);
end
fileOut = fullfile('SurrogateMaps',fileOut);
save(fileOut,'nullMaps','rho','d0');
fprintf(1,'Saved results out to %s\n',fileOut);

end
