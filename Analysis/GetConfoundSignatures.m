% GetConfoundSignatures
%-------------------------------------------------------------------------------
% Idea is to loop through a bunch of 'confounds' and determine their
% enrichment signatures -- can then be visualized alongside enrichment results
%-------------------------------------------------------------------------------

% Global settings (human or mouse):
whatSpecies = 'mouse';
resultsTables = struct();
doSave = true; % whether to save out to a .mat file at the end

%-------------------------------------------------------------------------------
% Define the parameter sets (to loop over later)
%-------------------------------------------------------------------------------
paramsLoop = struct();
switch whatSpecies
case 'human'
    paramsLoop.default = GiveMeDefaultParams(whatSpecies); % all defaults
    paramsLoop.noNorm = GiveMeDefaultParams(whatSpecies); % all defaults
    paramsLoop.noNorm.g.normalizationInternal = 'none';
    paramsLoop.normG = GiveMeDefaultParams(whatSpecies); % all defaults
    paramsLoop.normG.g.normalizationInternal = 'none';
    paramsLoop.normG.g.normalizationGene = 'zscore';
    paramsLoop.normGR = GiveMeDefaultParams(whatSpecies); % all defaults
    paramsLoop.normGR.g.normalizationInternal = 'none';
    paramsLoop.normGR.g.normalizationGene = 'zscore';
    paramsLoop.normGR.g.normalizationRegion = 'zscore';
case 'mouse'
    paramsLoop.default = GiveMeDefaultParams(whatSpecies); % all defaults
    paramsLoop.noNorm = paramsLoop.default;
    paramsLoop.noNorm.g.normalizationGene = 'none';
    paramsLoop.noNorm.g.normalizationRegion = 'none';
    paramsLoop.normG = paramsLoop.default;
    paramsLoop.normG.g.normalizationGene = 'zscore';
    paramsLoop.normG.g.normalizationRegion = 'none';
    paramsLoop.normR = paramsLoop.default;
    paramsLoop.normR.g.normalizationGene = 'none';
    paramsLoop.normR.g.normalizationRegion = 'zscore';
    paramsLoop.normGR = paramsLoop.default;
    paramsLoop.normGR.g.normalizationGene = 'zscore';
    paramsLoop.normGR.g.normalizationRegion = 'zscore';
end
paramSetLabels = fieldnames(paramsLoop);
numParamSets = length(paramSetLabels);
fprintf(1,'Comparing %u parameter sets for %s\n',numParamSets,whatSpecies);

%-------------------------------------------------------------------------------
% Human -- coexpression increases/decreases with distance
%-------------------------------------------------------------------------------
whatEdgeMeasure = 'distance';
onlyOnEdges = false;
correctDistance = false;
corrType = 'Spearman';
whatNull = 'randomGene';
numNulls = 100;

for i = 1:numParamSets
    paramsHere = paramsLoop.(paramSetLabels{i});
    absType = 'neg';
    tableName = sprintf('distance_neg_%s',paramSetLabels{i});
    [resultsTables.(tableName),~] = EdgeEnrichment(whatEdgeMeasure,onlyOnEdges,...
                    correctDistance,absType,corrType,whatNull,numNulls,whatSpecies,paramsHere);

    absType = 'pos';
    tableName = sprintf('distance_pos_%s',paramSetLabels{i});
    [resultsTables.(tableName),~] = EdgeEnrichment(whatEdgeMeasure,onlyOnEdges,...
                    correctDistance,absType,corrType,whatNull,numNulls,whatSpecies,paramsHere);
end

%===============================================================================
% Nodal correlations (in whole brain/cortex?)
%===============================================================================
switch whatSpecies
case 'mouse'
    structFilters = {'all','cortex'}; % compare both
case 'human'
    structFilters = {'cortex'}; % just look in cortex
end
numStructFilters = length(structFilters);
corrType = 'Spearman';

for k = 1:numStructFilters
    structFilterHere = structFilters{k};
    for i = 1:numParamSets
        paramsHere = paramsLoop.(paramSetLabels{i});

        % Highest variance:
        enrichWhat = 'varExpression';
        tableName = sprintf('varExpression_%s_%s',structFilterHere,paramSetLabels{i});
        [resultsTables.(tableName),~] = NodeSimpleEnrichment(enrichWhat,structFilterHere,corrType,whatSpecies,paramsHere);

        % Correlate with nodal degree:
        enrichWhat = 'degree';
        tableName = sprintf('degree_%s_%s',structFilterHere,paramSetLabels{i});
        [resultsTables.(tableName),~] = NodeSimpleEnrichment(enrichWhat,structFilterHere,corrType,whatSpecies,paramsHere);

        % Vary with dominant PC:
        enrichWhat = 'genePC';
        tableName = sprintf('genePC1_%s_%s',structFilterHere,paramSetLabels{i});
        [resultsTables.(tableName),~] = NodeSimpleEnrichment(enrichWhat,structFilterHere,corrType,whatSpecies,paramsHere);
    end
end

%===============================================================================
% Cortex/subcortex (mouse)
%===============================================================================
if strcmp(whatSpecies,'mouse')
    structFilter = 'all';
    corrType = 'Spearman';
    % Vary with dominant PC:
    enrichWhat = 'isocortex';

    for i = 1:numParamSets
        tableName = sprintf('cortex_diff_%s',paramSetLabels{i});
        [resultsTables.(tableName),~] = NodeSimpleEnrichment(enrichWhat,structFilter,corrType,whatSpecies,paramsHere);
    end
end

%===============================================================================
% Now for some saving
%===============================================================================
fileName = sprintf('resultsTables_%s.mat',whatSpecies);
fileName = fullfile('DataOutputs',fileName);
save(fileName,'resultsTables','whatSpecies');
fprintf(1,'Saved all results tables to ''%s''!!\n',fileName);

%===============================================================================
% Now for some plotting:
%===============================================================================
theThreshold = 0.2;

%===============================================================================
% Plot with human literature results:
%===============================================================================
resultsTablesLiterature = LiteratureLook(whatSpecies,theThreshold,false);

% Combine:
resultsTablesTogether = struct();
tableLabels = struct();
resultsNames = fieldnames(resultsTables);
numResultsHere = length(resultsNames);
litResultsNames = fieldnames(resultsTablesLiterature);
numLitResultsHere = length(litResultsNames);
% Results from our in-house enrichment:
for i = 1:numResultsHere
    resultsTablesTogether.(resultsNames{i}) = resultsTables.(resultsNames{i});
    tableLabels.(resultsNames{i}) = 1;
end
% Results from our literature survey:
for i = 1:numLitResultsHere
    resultsTablesTogether.(litResultsNames{i}) = resultsTablesLiterature.(litResultsNames{i});
    tableLabels.(litResultsNames{i}) = 2;
end

%-------------------------------------------------------------------------------
% Plot just the 'confound' signatures:
% PlotEnrichmentTables(resultsTables,theThreshold,whatSpecies);
% Plot a combined representation of 'confound' and 'data' signatures:
PlotEnrichmentTables(resultsTablesTogether,theThreshold,whatSpecies);
