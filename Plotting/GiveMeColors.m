function colorMatrix = GiveMeColors(whatColors)

switch whatColors
case 'twoGOCategories'
    colorMatrix = [161,222,240; 133,22,87; 127.5,127.5,127.5]/255;
    % colorMatrix = [136,233,154; 214,6,26]/255;
case 'mouseHuman'
    colorMatrix = [161,222,240; 17,93,82]/255;
case 'nullModels'
    colorMatrix = {[69,193,124]/255, [141,16,43]/255, [17,204,220]/255};
end

end
