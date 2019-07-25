function colorMatrix = GiveMeColors(whatColors)

switch whatColors
case 'twoGOCategories'
    colorMatrix = [161,222,240; 133,22,87]/255;
    % colorMatrix = [136,233,154; 214,6,26]/255;
case 'mouseHuman'
    colorMatrix = [161,222,240;17,93,82]/255;
end

end
