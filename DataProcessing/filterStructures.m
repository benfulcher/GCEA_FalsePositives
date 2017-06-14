function [A,geneData,structInfo,keepStruct] = filterStructures(structFilter,structInfo,A,geneData)

if isempty(structInfo)
    dataFile = '/Users/benfulcher/GoogleDrive/Work/CurrentProjects/CellTypesMouse/Code/AllenGeneDataset_19419.mat';
    load(dataFile,'structInfo');
end

switch structFilter
case 'cortex'
    keepStruct = strcmp(structInfo.divisionLabel,'Isocortex');
    geneData = geneData(keepStruct,:);
    structInfo = structInfo(keepStruct,:);
    A = A(keepStruct,keepStruct);
case 'all'
    keepStruct = 1:height(structInfo);
end

end
