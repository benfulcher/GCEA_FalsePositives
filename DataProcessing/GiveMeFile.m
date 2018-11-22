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
case 'HumanGene_cust100_normSeparate'
    % Cortex and subcortex normalized separately:
    dataFile = '100DS220scaledRobustSigmoidNSGDSQC1LcortexSubcortexSEPARATE_ROI_NOdistCorrSurfaceANDEuclidean.mat';
case 'HumanGene_cust100_normTogether'
    % Cortex and subcortex normalized together:
    dataFile = '100DS220scaledRobustSigmoidNSGDSQC1LcortexSubcortex_ROI_NOdistCorrSurfaceANDEuclidean.mat';
case 'HumanGene_HCP'
    % HCP parcellation of the cortex:
    dataFile = '100DS360scaledRobustSigmoidNSGDSQC1Lcortex_ROI_NOdistCorrSurface.mat';
otherwise
    error('Unknown file label %s',fileLabel);
end

end
