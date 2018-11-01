function dataFile = GiveMeFile(fileLabel)
% Returns a specific filename or path to external dependencies that vary across systems

switch fileLabel
case 'AllenMouseGene'
    if ismac
        dataFile = '/Users/benfulcher/DropboxSydneyUni/CompletedProjects/CellTypesMouse/Code/Data/AllenGeneDataset_19419.mat';
    else
        dataFile = 'MouseData/AllenGeneDataset_19419.mat';
    end
case 'EnrichmentToolbox'
    if ismac
        dataFile = '/Users/benfulcher/DropboxSydneyUni/CodeToolboxes/MatlabEnrichment/';
    else
        dataFile = '~/GeneEnrichment/';
    end
otherwise
    error('Unknown file label %s',fileLabel);
end

end
