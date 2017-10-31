function [A,geneData,structInfo,keepStruct] = filterStructures(structFilter,structInfo,A,geneData)

if nargin < 2 || isempty(structInfo)
    dataFile = '/Users/benfulcher/GoogleDrive/Work/CurrentProjects/CellTypesMouse/Code/Data/AllenGeneDataset_19419.mat';
    fprintf(1,'Importing structure information for mouse from file: %s\n',dataFile);
    load(dataFile,'structInfo');
end

if ismember('isCortex',structInfo.Properties.VariableNames)
    isHuman = true;
    fprintf(1,'Human\n');
else
    isHuman = false;
    fprintf(1,'Mouse\n');
end

%-------------------------------------------------------------------------------

switch structFilter
case {'isocortex','cortex'}
    if isHuman
        keepStruct = structInfo.isCortex;
    else
        keepStruct = strcmp(structInfo.divisionLabel,'Isocortex');
        end
    geneData = geneData(keepStruct,:);
    structInfo = structInfo(keepStruct,:);
    A = A(keepStruct,keepStruct);
case 'all'
    keepStruct = 1:height(structInfo);
end

end
