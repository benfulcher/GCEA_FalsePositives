function [f,ax] = PlotScatterPlot(structInfo,xData,yData,xLabel,yLabel,whatCorr,newFigure)
%-------------------------------------------------------------------------------
% In this case you provide values rather than just an ordering
%-------------------------------------------------------------------------------

if nargin < 6
    whatCorr = 'Pearson';
end
if nargin < 7
    newFigure = 1;
end

numRegions = height(structInfo);
if numRegions < 100
    labelCortex = 1;
else
    labelCortex = 0;
end

%-------------------------------------------------------------------------------
if newFigure
    f = figure('color','w');
else
    f = gcf;
end
hold on; ax = gca;

%-------------------------------------------------------------------------------
% Also label regions by their broad cortical area:
if labelCortex
    [areaLabels,labelInd,labelNames] = LabelCorticalAreas(structInfo.acronym);

    areaColors = BF_getcmap('pastel2',max(labelInd),1);
    plotHandles = cell(max(labelInd),1);
    for i = 1:max(labelInd)
        plotHandles{i} = plot(xData(labelInd==i),yData(labelInd==i),...
                        'o','MarkerFaceColor',areaColors{i},'MarkerSize',12);
    end

    % Add labels:
    xDataRange = range(xData);
    for i = 1:numRegions
        text(xData(i)+0.04*xDataRange,yData(i),structInfo.acronym{i},...
                            'color',brighten(areaColors{labelInd(i)},-0.7))
    end

    % Add a legend
    legend([plotHandles{:}],labelNames,'Location','best')
else
    % Scatter plot with the region colors:
    dotColors = arrayfun(@(x)rgbconv(structInfo.color_hex_triplet{x})',...
                                            1:numRegions,'UniformOutput',0);
    dotColors = [dotColors{:}]';

    nodeSize = 50;
    scatter(xData,yData,nodeSize,dotColors,'fill','MarkerEdgeColor','k')

    % Add labels:
    xDataRange = range(xData);
    for i = 1:numRegions
        text(xData(i)+0.04*xDataRange,yData(i),structInfo.acronym{i},...
                            'color',brighten(dotColors(i,:),-0.3))
    end

    % And for major divisions:
    divisionLabels = categorical(structInfo.divisionLabel);
    theDivisions = unique(divisionLabels);
    numDivisions = length(theDivisions);
    for i = 1:numDivisions
        % Put each major region in the center of those points
        centrePoint = [mean(xData(divisionLabels==theDivisions(i))),mean(yData(divisionLabels==theDivisions(i)))];
        find_1 = find(divisionLabels==theDivisions(i),1);
        text(centrePoint(1),centrePoint(2),char(theDivisions(i)), ...
                    'color','k','FontSize',14,'BackgroundColor',dotColors(find_1,:))
    end
end

%-------------------------------------------------------------------------------
% Cosmetics:
xlabel(xLabel,'interpreter','none')
ylabel(yLabel,'interpreter','none')

% Add title:
isGood = ~isnan(yData) & ~isnan(xData);
[rho,pVal] = corr(xData(isGood),yData(isGood),'type',whatCorr);
title(sprintf('%s r = %.3f (p = %.3g)',whatCorr,rho,pVal))

end
