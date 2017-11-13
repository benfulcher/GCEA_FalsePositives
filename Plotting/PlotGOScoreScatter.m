function PlotGOScoreScatter(resultsTab1,resultsTab2)
% Plots correlation in mean scores for all GO categories, from two
% different analyses
%-------------------------------------------------------------------------------

[~,ia,ib] = intersect(resultsTab1.GOID,resultsTab2.GOID);
xData = resultsTab1.meanScore(ia);
yData = resultsTab2.meanScore(ib);

f = figure('color','w');
plot(xData,yData,'.k')
r = corr(resultsTab1.meanScore(ia),resultsTab2.meanScore(ib),'rows','pairwise');
title(sprintf('%u categories, r = %.3f',length(ia),r))
axis square
f.Position = [1000,1027,432,311];

end
