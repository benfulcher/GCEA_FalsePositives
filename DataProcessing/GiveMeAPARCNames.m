function regionStruct = GiveMeAPARCNames()
% Loads in Aurina's data on APARC names, and outputs

fid = fopen('ROInames_aparcasegBen.txt');
S = textscan(fid,'%u%s');
fclose(fid);
ID = S{1};
Label = S{2};

% Add metadata:
isRightCortex = cellfun(@(x)strcmp(x(1:6),'ctx-rh'),Label);
isRightSubcortex = cellfun(@(x)strcmp(x(1:6),'Right-'),Label);
isLeftCortex = cellfun(@(x)strcmp(x(1:6),'ctx-lh'),Label);
isLeftSubcortex = cellfun(@(x)strcmp(x(1:5),'Left-'),Label);

isRight = (isRightCortex | isRightSubcortex);
isLeft = (isLeftCortex | isLeftSubcortex);
isCortex = (isRightCortex | isLeftCortex);
isSubcortex = (isRightSubcortex | isLeftSubcortex);

regionStruct = table(ID,Label,isLeft,isRight,isCortex,isSubcortex);

end
