function fileSimpleEnrichmentResults = GiveMeSimpleEnrichmentOutputFile(params,whatProperty)

% Construct a basic file name:
fileSimpleEnrichmentResults = sprintf('SimpleEnrichment_%s_%u_%s-%s.mat',...
        whatProperty,params.e.numNullSamples,params.humanOrMouse,params.structFilter);

% Put it in a separate directory of ensemble-based nulls:
fileSimpleEnrichmentResults = fullfile('DataOutputs',fileSimpleEnrichmentResults);

end
