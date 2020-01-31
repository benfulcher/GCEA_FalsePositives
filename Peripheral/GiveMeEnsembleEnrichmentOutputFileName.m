function fileNullEnsembleResults = GiveMeEnsembleEnrichmentOutputFileName(params,doOld)
% Get the filename of interest
%-------------------------------------------------------------------------------
% Don't use the old results by default!
if nargin < 2
    doOld = false;
end
%-------------------------------------------------------------------------------
% Need to find precomputed null results (from running ComputeAllCategoryNulls):
if doOld
    warning('Using old data! Please assure me this is only temporary...')
    fileNullEnsembleResults = sprintf('RandomNull_%u_%s-%s_%s_%s_%s.mat',...
                params.e.numNullSamples,params.humanOrMouse,...
                params.structFilter,params.e.whatEnsemble,...
                params.e.whatCorr,params.e.aggregateHow);
else
    fileNullEnsembleResults = sprintf('PhenotypeNulls_%u_%s-%s_%s_%s_%s.mat',...
                params.e.numNullSamples,params.humanOrMouse,...
                params.structFilter,params.e.whatEnsemble,...
                params.e.whatCorr,params.e.aggregateHow);
end

end
