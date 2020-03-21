function NullEnrichmentTogether(whatSpecies,doLog)
% Investigate the distribution of category-level false-positive (corrected) rates
% across different ensembles of null phenotypes.
%-------------------------------------------------------------------------------
if nargin < 1
    whatSpecies = 'human';
end
params = GiveMeDefaultParams(whatSpecies);
if nargin < 2
    doLog = true; % logarithmic vertical axis
end
%-------------------------------------------------------------------------------
doDisplay = false; % info about each CFPR result

nullGOTables = struct();

% (i) REFERENCE: random map ensemble with independent spatial shuffling
params.g.whatSurrogate = 'randomMap';
params.nulls.customShuffle = 'independentSpatialShuffle';
nullGOTables.reference = SurrogateEnrichmentProcess(params,doDisplay);

% (ii) SBP-random: random map ensemble
params.nulls.customShuffle = 'none';
nullGOTables.SBPrandom = SurrogateEnrichmentProcess(params,doDisplay);

% (iii) SBP-spatial: ensemble of spatially autocorrelated maps
params.g.whatSurrogate = 'spatialLag';
nullGOTables.SBPspatial = SurrogateEnrichmentProcess(params,doDisplay);

% nullGOTables.coordSpatialRandom = SurrogateEnrichmentProcess(whatSpecies,params.nulls.numNullsCFPR,'randomUniform','coordinatedSpatialShuffle');

%-------------------------------------------------------------------------------
% Extract sum under significant values:
sumUnderSigValues = structfun(@(x)x.sumUnderSig/params.nulls.numNullsCFPR,...
                                            nullGOTables,'UniformOutput',false);

%-------------------------------------------------------------------------------
% Distribution of sumUnderSig:
fields = fieldnames(nullGOTables);
sumUnderSigCell = struct2cell(sumUnderSigValues);

permutationKeep = [1,2,3];

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
BF_JitteredParallelScatter(sumUnderSigCell,false,true,true,extraParams);

ax = gca();
ax.XTick = 1:3;
ax.XTickLabel = fields(permutationKeep);
ylabel('Category false-positive rate (CFPR)')
% title('Distribution over GO categories')
f = gcf();
f.Position = [1000        1121         273         217];

% Save:
fileName = fullfile('OutputPlots',sprintf('CFPR_distributions_%s.svg',whatSpecies));
saveas(f,fileName,'svg')
fprintf(1,'Saved to %s\n',fileName);

end
