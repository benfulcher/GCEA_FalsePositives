function yDataCorrected = BF_PlotQuantiles(xData,yData,numBins,alsoScatter,makeNewFigure,customColor,doStds)
% Plots x-y scatter, but with mean of y plotted in quantiles of x
% Outputs yData, corrected for the quantile means
% Ben Fulcher
%-------------------------------------------------------------------------------

if nargin < 3 || isempty(numBins)
    numBins = 10;
end
numThresholds = numBins + 1;
if nargin < 4
    alsoScatter = false;
end
if nargin < 5
    makeNewFigure = false;
end
if nargin < 6
    customColor = 'r';
end
if nargin < 7
    doStds = true;
end

%-------------------------------------------------------------------------------
% Filter out NaNs:
goodBoth = (isfinite(xData) & isfinite(yData));
if ~any(goodBoth)
    error('No good data');
elseif any(~goodBoth)
    xData = xData(goodBoth);
    yData = yData(goodBoth);
    fprintf(1,'Removed %u bad samples from x/y data\n',sum(~goodBoth));
end

xThresholds = arrayfun(@(x)quantile(xData,x),linspace(0,1,numThresholds));
xThresholds(end) = xThresholds(end) + eps; % make sure all data included in final bin
yMeans = arrayfun(@(x)mean(yData(xData>=xThresholds(x) & xData < xThresholds(x+1))),1:numThresholds-1);
yStds = arrayfun(@(x)std(yData(xData>=xThresholds(x) & xData < xThresholds(x+1))),1:numThresholds-1);

% ------------------------------------------------------------------------------
% Plot:
if makeNewFigure
    f = figure('color','w'); box('on');
end
hold('on')
theColor = customColor;
theStyle = '-';
theLineWidth = 2;

if alsoScatter
    plot(xData,yData,'.k');
end

for k = 1:numThresholds-1
    plot(xThresholds(k:k+1),ones(2,1)*yMeans(k),'LineStyle',theStyle,'LineWidth',theLineWidth,'Color',theColor)
    plot(mean(xThresholds(k:k+1)),yMeans(k),'o','MarkerSize',5,'LineStyle',theStyle,'LineWidth',theLineWidth,'Color',theColor)
    if doStds
        plot(xThresholds(k:k+1),ones(2,1)*(yMeans(k)+yStds(k)),'LineStyle','--','LineWidth',theLineWidth,'Color',theColor)
        plot(xThresholds(k:k+1),ones(2,1)*(yMeans(k)-yStds(k)),'LineStyle','--','LineWidth',theLineWidth,'Color',theColor)
    end
end

%-------------------------------------------------------------------------------
% Correct the yData
yDataCorrected = nan(size(yData));
for p = 1:numThresholds-1
    inBin = (xData>=xThresholds(p) & xData < xThresholds(p+1));
    yDataCorrected(inBin) = yData(inBin) - yMeans(p);
end


end
