function PlotColorMatrix(geneData,structInfo,colorLabelsWhere,myColorMap)
% Plots a colored data matrix, with the mouse connectome regions labeled

% ------------------------------------------------------------------------------
% Check inputs:
if nargin < 2
    params = GiveMeDefaultParams('mouse');
    [geneData,geneInfo,structInfo] = LoadMeG(params.g);
end
if nargin < 3 || isempty(colorLabelsWhere)
    colorLabelsWhere = 'left'; % left, bottom, both
end
if nargin < 5 || isempty(myColorMap)
    myColorMap = [flipud(BF_getcmap('blues',9));1,1,1;BF_getcmap('reds',9)];
end

%-------------------------------------------------------------------------------
% Additional custom parameters:
rectThickness = arrayfun(@(x)size(geneData,x)/50,1:2);
if length(rectThickness)==1
    rectThickness = ones(2,1)*rectThickness;
end
labelInd = false; % Labels individual regions rather than major regions
plotBoundaries = [false,false];

% ------------------------------------------------------------------------------
% Plot the data matrix
f = figure('color','w');
ax = gca();
hold('on');
BF_imagesc(geneData);
colormap(myColorMap)

%-------------------------------------------------------------------------------
% Y-axis labels in the middle of each contiguous region of major region labels
majorLabels = structInfo.divisionLabel;
[~,ia,ib] = unique(majorLabels,'stable');
ia_mid = [ia;length(majorLabels)];
ia_mid = floor(mean([ia_mid(1:end-1),ia_mid(2:end)],2));
ax.YTick = ia_mid;
ax.YTickLabel = majorLabels(ia_mid);

%-------------------------------------------------------------------------------
% Add rectangles labeling major brain regions, and individual colors
% ------------------------------------------------------------------------------
if ismember('color_hex_triplet',structInfo.Properties.VariableNames)
    for j = 1:height(structInfo)
        colorHere = rgbconv(structInfo.color_hex_triplet{j});
        rectangle('Position',[1-rectThickness(2),j,rectThickness(2),1], ...
                    'FaceColor',colorHere,'EdgeColor',colorHere)
    end
end

% ------------------------------------------------------------------------------
% Add separator black lines:
% ------------------------------------------------------------------------------
lineWidth = 1.2;
if ismember(colorLabelsWhere,{'both','left'})
    for j = 1:length(ia)
        plot([1-rectThickness(2),1],ones(2,1)*ia(j),'k','LineWidth',lineWidth)
    end
    % Bottom one:
    plot([1-rectThickness(2),1],ones(2,1)*1,'k','LineWidth',lineWidth)
    % Top one:
    plot([1-rectThickness(2),1],ones(2,1)*size(geneData,1)+1,'k','LineWidth',lineWidth)
    % Left one:
    plot((1-rectThickness(2))*ones(2,1),[1,size(geneData,1)+1],'k','LineWidth',lineWidth)
    % Right one:
    plot(ones(2,1),[1,size(geneData,1)+1],'k','LineWidth',lineWidth)
end
if ismember(colorLabelsWhere,{'both','right'})
    for j = 1:length(ia)
        plot(size(geneData,2)+1+[0,rectThickness(2)],ones(2,1)*ia(j),'k','LineWidth',lineWidth)
    end
    % Top one:
    plot(size(geneData,2)+1+[0,rectThickness(2)],ones(2,1)*size(geneData,1)+1,'k','LineWidth',lineWidth)
    % Left one:
    plot((size(geneData,2)+1)*ones(2,1),[1,size(geneData,1)+1],'k','LineWidth',lineWidth)
    % Right one:
    plot((size(geneData,2)+1+rectThickness(2))*ones(2,1),[1,size(geneData,1)+1],'k','LineWidth',lineWidth)
    % Bottom one:
    plot(size(geneData,2)+1+[0,rectThickness(2)],ones(2,1),'k','LineWidth',lineWidth)
end
if ismember(colorLabelsWhere,{'both','bottom'})
    for j = 1:length(ia)
        plot(ones(2,1)*(size(geneData,1)-ia(j)+2),[1-rectThickness(1),1],'k','LineWidth',lineWidth)
    end
    % Leftmost one:
    plot(ones(2,1),[1-rectThickness(1),1],'k','LineWidth',lineWidth)
    % Bottom one:
    plot([1,size(geneData,2)+1],(1-rectThickness(1))*ones(2,1),'k','LineWidth',lineWidth)
    % Right one:
    plot(ones(2,1)*(size(geneData,2)+1),[1-rectThickness(1),1],'k','LineWidth',lineWidth)
    % Top one:
    plot([1,size(geneData,2)+1],ones(2,1),'k','LineWidth',lineWidth)
end
if ismember(colorLabelsWhere,{'both','top'})
    if ~(labelInd || isPermuted)
        for j = 1:length(ia)
            plot(ones(2,1)*(size(geneData,2)-ia(j)+2),size(geneData,1)+1+[0,rectThickness(1)], ...
                            'k','LineWidth',lineWidth)
        end
    end
    % Leftmost one:
    plot(ones(2,1),size(geneData,1)+1+[0,rectThickness(1)],'k','LineWidth',lineWidth)
    % Rightmost one:
    plot(ones(2,1)*(size(geneData,2)+1),size(geneData,1)+1+[0,rectThickness(1)],'k','LineWidth',lineWidth)
    % Bottom one:
    plot([1,size(geneData,2)+1],(size(geneData,1)+1)*ones(2,1),'k','LineWidth',lineWidth)
    % Top one:
    plot([1,size(geneData,2)+1],(size(geneData,1)+1+rectThickness(1))*ones(2,1),'k','LineWidth',lineWidth)
end

% ------------------------------------------------------------------------------
% Adjust axes to see the labeling:
% ------------------------------------------------------------------------------
scaleFactor = 0.08; % to see the little bit extra to capture the line thickness
% First set with the scale factor
switch colorLabelsWhere
case 'both'
    ax.XLim = [-rectThickness(2)*(scaleFactor),size(geneData,2)+rectThickness(2)*scaleFactor];
    ax.YLim = [-rectThickness(1)*(scaleFactor),size(geneData,1)+rectThickness(1)*scaleFactor];
case 'left'
    ax.XLim = [-rectThickness(2)*(scaleFactor),size(geneData,2)+rectThickness(2)*scaleFactor];
end
if ismember(colorLabelsWhere,{'both','left'})
    ax.XLim = [1-rectThickness(2)*(1+scaleFactor),size(geneData,2)+1];
end
if ismember(colorLabelsWhere,{'both','right'})
    ax.XLim(2) = size(geneData,2) + 2 + rectThickness(2)*(1+scaleFactor);
end
if ismember(colorLabelsWhere,{'both','bottom'})
    ax.YLim(1) = -rectThickness(1)*(1+scaleFactor);
end
if ismember(colorLabelsWhere,{'both','top'})
    ax.YLim(2) = size(geneData,1) + 2 + rectThickness(1)*(1+scaleFactor);
end

% Remove tick marks:
set(gca,'TickLength',[0,0])

end
