function NullEnrichmentTogether(whatSpecies,doLog)

if nargin < 1
    whatSpecies = 'human';
end
if nargin < 2
    doLog = true;
end

%-------------------------------------------------------------------------------
numNullSamples_surrogate = 10000; % (SurrogateGOTables_10000_*.mat)

%-------------------------------------------------------------------------------

nullGOTables = struct();
nullGOTables.spatialRandom = SurrogateEnrichmentProcess(whatSpecies,numNullSamples_surrogate,'randomUniform','');
% nullGOTables.coordSpatialRandom = SurrogateEnrichmentProcess(whatSpecies,numNullSamples_surrogate,'randomUniform','coordinatedSpatialShuffle');
nullGOTables.indSpatialRandom = SurrogateEnrichmentProcess(whatSpecies,numNullSamples_surrogate,'randomUniform','independentSpatialShuffle');
nullGOTables.spatialLag = SurrogateEnrichmentProcess(whatSpecies,numNullSamples_surrogate,'spatialLag','');

%-------------------------------------------------------------------------------
% Extract sum under significant values:
sumUnderSigValues = structfun(@(x)x.sumUnderSig/numNullSamples_surrogate,nullGOTables,'UniformOutput',false);

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
