% Quick look at whether reciprocal connections have correlated edge values
%-------------------------------------------------------------------------------

connectomeSource = 'Oh';
pThreshold = 0.05;
whatHemispheres = 'right';
justCortex = false;
whatEdgeMeasure = 'bin_edgeBet';
onlyOnEdges = true; % whether to put values only on existing edges
                    % (rather than all node pairs for some measures)

[A_bin,regionStruct,adjPVals] = GiveMeAdj(connectomeSource,pThreshold,true,whatHemispheres,justCortex);
edgeMetric = edge_betweenness_bin(A_bin);
% ktot = sum(A_bin,1)' + sum(A_bin,2);
% product = ktot*ktot';
% product(A_bin == 0) = 0;
% edgeMetric = product;

%-------------------------------------------------------------------------------
isRecip = (A_bin==1 & A_bin'==1);
[i,j] = find(isRecip);
isLower = (i<j);
isUpper = (i>j);
fisLower = find(isLower);
numPairs = sum(isLower);
matchPairs = zeros(numPairs,2);
for p = 1:numPairs
    matchPairs(p,:) = [edgeMetric(i(fisLower(p)),j(fisLower(p))),edgeMetric(j(fisLower(p)),i(fisLower(p)))];
end
f = figure('color','w');
plot(matchPairs(:,1),matchPairs(:,2),'.k')
r = corr(matchPairs(:,1),matchPairs(:,2));
title(r)
