function dataFile = GiveMeFile(fileLabel)
% Returns a specific filename or path to external dependencies that vary across systems
%-------------------------------------------------------------------------------

switch fileLabel
case 'AllenMouseGene'
    dataFile = fullfile('MouseData','AllenGeneDataset_19419.mat');
case 'EnrichmentToolbox'
    if ismac
        dataFile = '~/DropboxSydneyUni/CodeToolboxes/GeneCategoryEnrichment/CodeRepo/';
    else
        % Can install by git clone git@GeneSetEnrichmentAnalysis.git
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
