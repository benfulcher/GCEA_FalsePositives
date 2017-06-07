function RegionPlot(structInfo,structData)

numNodes = height(structInfo);

% ------------------------------------------------------------------------------
% Plot barcodes
% ------------------------------------------------------------------------------
f = figure('color','w'); box('on'); ax = gca;
for k = 1:numNodes
    colorHere = rgbconv(structInfo.color_hex_triplet{k});
    rectangle('Position',[0,numNodes-(k-0.5),structData(k),1], ...
                'FaceColor',colorHere,'EdgeColor',colorHere)
end
xlabel('degree')
ax.YTick = 1:numNodes;
ax.YLim = 0.5+[0,numNodes];
ax.YTickLabel = structInfo.acronym;

end
