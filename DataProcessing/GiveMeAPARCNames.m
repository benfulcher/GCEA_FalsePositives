function regionStruct = GiveMeAPARCNames()
% Loads in Aurina's data on APARC names, and outputs

fid = fopen('ROInames_aparcasegBen.txt');
S = textscan(fid,'%u%s');
fclose(fid);
ID = S{1};
acronym = S{2};

% Add metadata:
isRightCortex = cellfun(@(x)strcmp(x(1:6),'ctx-rh'),acronym);
isRightSubcortex = cellfun(@(x)strcmp(x(1:6),'Right-'),acronym);
isLeftCortex = cellfun(@(x)strcmp(x(1:6),'ctx-lh'),acronym);
isLeftSubcortex = cellfun(@(x)strcmp(x(1:5),'Left-'),acronym);

isRight = (isRightCortex | isRightSubcortex);
isLeft = (isLeftCortex | isLeftSubcortex);
isCortex = (isRightCortex | isLeftCortex);
isSubcortex = (isRightSubcortex | isLeftSubcortex);

regionStruct = table(ID,acronym,isLeft,isRight,isCortex,isSubcortex);

end
