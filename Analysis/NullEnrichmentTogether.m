function NullEnrichmentTogether(whatSpecies,doLog)

if nargin < 1
    whatSpecies = 'mouse';
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
sumUnderSigCell = sumUnderSigCell([3,1,2])
% f = figure('color','w'); hold('on')
% histogram(sumUnderSigCell{3},'normalization','probability')
% histogram(sumUnderSigCell{1},'normalization','probability')
% histogram(sumUnderSigCell{2},'normalization','probability')

%-------------------------------------------------------------------------------
% Figure

if doLog
    % ***As log***:
    log10Data = cellfun(@(x)log10(x),sumUnderSigCell,'UniformOutput',false);
    log10Data = cellfun(@(x)x(isfinite(x)),log10Data,'UniformOutput',false);
    BF_JitteredParallelScatter(log10Data);
else
    % ***As linear***:
    BF_JitteredParallelScatter(sumUnderSigCell);
end

ax = gca();
ax.XTick = 1:3;
ax.XTickLabel = fields([3,1,2]);
ylabel('RPS (false positive rate)')
title('Distribution over GO categories')
f = gcf();
f.Position = [1000        1121         273         217];

end
