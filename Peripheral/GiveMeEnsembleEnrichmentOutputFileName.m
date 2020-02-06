function fileNullEnsembleResults = GiveMeEnsembleEnrichmentOutputFileName(params)
% Get the filename of interest
%-------------------------------------------------------------------------------
% Don't use the old results by default!
if nargin < 2
    doOld = false;
end
%-------------------------------------------------------------------------------
% Need to find precomputed null results (from running ComputeAllCategoryNulls):

% Try new:
fileNullEnsembleResults = sprintf('PhenotypeNulls_%u_%s-%s_%s_%s_%s.mat',...
            params.e.numNullSamples,params.humanOrMouse,...
            params.structFilter,params.e.whatEnsemble,...
            params.e.whatCorr,params.e.aggregateHow);

% Put it in a separate directory of ensemble-based nulls:
fileNullEnsembleResults = fullfile('EnsembleBasedNulls',fileNullEnsembleResults);

%-------------------------------------------------------------------------------
% What if it's not computed yet? See if there's an old one we can use...
%-------------------------------------------------------------------------------
if ~exist(fileNullEnsembleResults)
    warning('Using old data! Please assure me this is only temporary...')
    fileNullEnsembleResults = sprintf('RandomNull_%u_%s-%s_%s_%s_%s.mat',...
                params.e.numNullSamples,params.humanOrMouse,...
                params.structFilter,params.e.whatEnsemble,...
                params.e.whatCorr,params.e.aggregateHow);
    fileNullEnsembleResults = fullfile('EnsembleBasedNulls',fileNullEnsembleResults);
end

end
