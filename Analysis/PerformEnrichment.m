function GOTablePhenotype = PerformEnrichment(myPhenotype,whatSpecies,whatNullModel)
% PerformEnrichment  Compute enrichment in different GO categories according to
%                       a given null model
% Assumes the null samples have been precomputed using ComputeAllCategoryNulls

%---INPUTS:
% whatSpecies: 'mouse' or 'human'
% myPhenotype: the phenotype to assess for enrichment
% whatNullModel: the null model to use for enrichment

%-------------------------------------------------------------------------------
% INPUTS:
%-------------------------------------------------------------------------------
if nargin < 1
    myPhenotype = 'degree';
end
if nargin < 2
    whatSpecies = 'mouse';
end
params = GiveMeDefaultParams(whatSpecies);
if nargin < 3
    whatNullModel = 'randomMap'; % 'spatialLag'
end

whatCorr = 'Spearman';
aggregateHow = 'mean';

%-------------------------------------------------------------------------------
% Load null distributions using default settings:
%-------------------------------------------------------------------------------
% Check for precomputed null results (ComputeAllCategoryNulls):
numNullSamples = 20000;
fileNameDesired = sprintf('RandomNull_%u_%s-%s_%s_%s_%s.mat',numNullSamples,...
                            whatSpecies,params.g.structFilter,whatNullModel,whatCorr,aggregateHow);
nullDistributions = load(fileNameDesired,'GOTable');
GOTableNull = nullDistributions.GOTable;
% (The categoryScores variable is the distribution of null samples for each GO category)

%-------------------------------------------------------------------------------
% Compute the phenotype of interest
%-------------------------------------------------------------------------------
switch myPhenotype
case 'degree'
    myPhenotype = ComputeDegree(whatSpecies,'binary');
end

%-------------------------------------------------------------------------------
% Assign equivalent scores to phenotype
%-------------------------------------------------------------------------------
GOTablePhenotype = ComputeAllCategoryNulls(params,1,myPhenotype,whatCorr,aggregateHow,false);

% Check that we have the same GO category IDs in both cases:
if ~all(GOTableNull.GOID==GOTablePhenotype.GOID)
    error('Error matching GO Categories to precomputed null data...');
end
numCategories = height(GOTablePhenotype);

%-------------------------------------------------------------------------------
% Estimate p-values:
%-------------------------------------------------------------------------------
GOTablePhenotype = EstimatePVals(GOTableNull.categoryScores,[GOTablePhenotype.categoryScores{:}],'right',GOTablePhenotype);


numSig = sum(GOTablePhenotype.pValZCorr < params.e.sigThresh);
fprintf(1,'%u significant categories at p_corr < %.2f\n',numSig,params.e.sigThresh);

[geneData,geneInfo,structInfo] = LoadMeG(params.g);
ix_GO = ListCategories(geneInfo,GOTablePhenotype,20,'pValZ');

end
