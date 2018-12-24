function PlotGOScoreScatter(resultsTab1,resultsTab2,customField)
% Plots correlation in mean scores for all GO categories, from two
% different analyses
%-------------------------------------------------------------------------------

% Check inputs, set defaults:
if nargin < 3
    customField = {'meanScore','meanScore'};
end
%-------------------------------------------------------------------------------

% Make sure we investigate GO categories that are common to both analyses:
[~,ia,ib] = intersect(resultsTab1.GOID,resultsTab2.GOID);
xData = resultsTab1.(customField{1})(ia);
yData = resultsTab2.(customField{2})(ib);

% Compute correlation:
r = corr(xData,yData,'type','Spearman','rows','pairwise');

% Plot the scatter:
f = figure('color','w');
plot(xData,yData,'.k')
title(sprintf('%u categories, r = %.3f',length(ia),r))
axis square
f.Position = [1000,1027,432,311];

end
