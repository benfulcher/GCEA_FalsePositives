function NullComputation(params)
% Wrapper for running ComputeAllCategoryNulls using appropriate gene-expression
% data
%-------------------------------------------------------------------------------

if nargin < 1
    params = GiveMeDefaultParams();
end
%-------------------------------------------------------------------------------

% Double check the output filename is appropriate for the parameter settings used:
params.e.fileNameOut = GiveMeEnsembleEnrichmentOutputFileName(params);

% Get gene-expression data to feed in:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
geneDataStruct = struct();
geneDataStruct.expressionMatrix = geneData;
geneDataStruct.entrezIDs = geneInfo.entrez_id;

% Run the computation using code from the Enrichment Repository:
ComputeAllCategoryNulls(geneDataStruct,params.e,[],true,true);

end
