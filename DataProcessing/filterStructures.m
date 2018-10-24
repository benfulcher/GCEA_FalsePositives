function [A,geneData,structInfo,keepStruct] = filterStructures(structFilter,structInfo,A,geneData)
% filterStructures   Reduces data to a specified subset of structures

%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    structFilter = 'cortex';
    fprintf(1,'Filtering to include only cortical areas by default\n');
end
if nargin < 2 || isempty(structInfo)
    dataFile = '/Users/benfulcher/GoogleDrive/Work/CompletedProjects/CellTypesMouse/Code/Data/AllenGeneDataset_19419.mat';
    fprintf(1,'Importing structure information for mouse from file: %s\n',dataFile);
    load(dataFile,'structInfo');
end
if nargin < 3
    A = [];
end
if ismember('isCortex',structInfo.Properties.VariableNames)
    isHuman = true;
    fprintf(1,'Human\n');
else
    isHuman = false;
    fprintf(1,'Mouse\n');
end

%-------------------------------------------------------------------------------
if isHuman
    isCortex = structInfo.isCortex;
else % mouse
    isCortex = strcmp(structInfo.divisionLabel,'Isocortex');
end

switch structFilter
case {'isocortex','cortex'}
    keepStruct = isCortex;
    fprintf(1,'Filtering to consider %u cortical areas\n',sum(keepStruct));
case 'notCortex'
    keepStruct = ~isCortex;
    fprintf(1,'Filtering to consider %u non-cortical areas\n',sum(keepStruct));
case 'all'
    keepStruct = true(height(structInfo),1);
    fprintf(1,'Keeping all %u areas\n',sum(keepStruct));
end

%-------------------------------------------------------------------------------
% Do the filtering:
if ~all(keepStruct)
    geneData = geneData(keepStruct,:);
    structInfo = structInfo(keepStruct,:);
    if ~isempty(A)
        A = A(keepStruct,keepStruct);
    end
end

end
