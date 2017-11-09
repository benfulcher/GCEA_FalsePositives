function structInfo = GiveMeAPARCNames()
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
if ~all(~isRight==isLeft)
    error('Error assigning left/right');
end
isCortex = (isRightCortex | isLeftCortex);
isSubcortex = (isRightSubcortex | isLeftSubcortex);
if ~all(~isCortex==isSubcortex)
    error('Error assigning cortex/subcortex');
end

structInfo = table(ID,acronym,isLeft,isCortex);

end
