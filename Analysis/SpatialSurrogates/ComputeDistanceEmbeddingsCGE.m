%-------------------------------------------------------------------------------
% Get exponential length-scales of gene-expression similarity?
%-------------------------------------------------------------------------------
mouseOrHuman = 'mouse';
structFilter = 'all';

params = GiveMeDefaultParams(mouseOrHuman,structFilter);
distMat = GiveMeDistanceMatrix(params);
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
G = corr(geneData','rows','pairwise');
upperMask = triu(true(size(G)),+1);
f = figure('color','w');
xData = distMat(upperMask);
yData = G(upperMask);
[f_handle,Stats,c] = GiveMeFit(xData,yData,'exp0',false);
plot(distMat(upperMask),G(upperMask),'.k')
