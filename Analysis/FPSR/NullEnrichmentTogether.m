function NullEnrichmentTogether(whatSpecies,numNullSamplesSurrogate,doLog)
% Investigate the distribution of false-positive (corrected) enrichment across
% different null phenotype ensembles.
%-------------------------------------------------------------------------------
if nargin < 1
    whatSpecies = 'human';
end
params = GiveMeDefaultParams(whatSpecies);
if nargin < 2 || isempty(numNullSamplesSurrogate)
    numNullSamplesSurrogate = params.nulls.numNullsFPSR; % number of null maps to test against
end
% (SurrogateGOTables_10000_*.mat)
if nargin < 3
    doLog = true;
end
%-------------------------------------------------------------------------------

nullGOTables = struct();
nullGOTables.SBPrandom = SurrogateEnrichmentProcess(whatSpecies,numNullSamplesSurrogate,'randomUniform','');
% nullGOTables.coordSpatialRandom = SurrogateEnrichmentProcess(whatSpecies,numNullSamplesSurrogate,'randomUniform','coordinatedSpatialShuffle');
nullGOTables.reference = SurrogateEnrichmentProcess(whatSpecies,numNullSamplesSurrogate,'randomUniform','independentSpatialShuffle');
nullGOTables.SBPspatial = SurrogateEnrichmentProcess(whatSpecies,numNullSamplesSurrogate,'spatialLag','');

%-------------------------------------------------------------------------------
% Extract sum under significant values:
sumUnderSigValues = structfun(@(x)x.sumUnderSig/numNullSamplesSurrogate,nullGOTables,'UniformOutput',false);

%-------------------------------------------------------------------------------
% Distribution of sumUnderSig:
fields = fieldnames(nullGOTables);
sumUnderSigCell = struct2cell(sumUnderSigValues);

permutationKeep = [2,1,3];

sumUnderSigCell = sumUnderSigCell(permutationKeep);

% f = figure('color','w'); hold('on')
% histogram(sumUnderSigCell{3},'normalization','probability')
% histogram(sumUnderSigCell{1},'normalization','probability')
% histogram(sumUnderSigCell{2},'normalization','probability')

%-------------------------------------------------------------------------------
% Figure
extraParams = struct();
extraParams.doLog = true;
extraParams.customSpot = '';
extraParams.theColors = GiveMeColors('nullModels');
BF_JitteredParallelScatter(sumUnderSigCell,false,true,true,extraParams)

ax = gca();
ax.XTick = 1:3;
ax.XTickLabel = fields(permutationKeep);
ylabel('False-positive significance rate (FPSR)')
title('Distribution over GO categories')
f = gcf();
f.Position = [1000        1121         273         217];

end
