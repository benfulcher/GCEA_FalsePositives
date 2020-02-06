function RunConventionalEnrichment()
% Runs all simple enrichment analyses for mouse/human cortex/whole-brain
%-------------------------------------------------------------------------------

speciesStructure = {'human','cortex';
                    'mouse','all';
                    'mouse','cortex'};
numAnalyses = size(speciesStructure,1);

for i = 1:numAnalyses
    params = GiveMeDefaultParams(speciesStructure{i,1},speciesStructure{i,2});
    NodeSimpleEnrichment(params,'degree',true);
end

end
