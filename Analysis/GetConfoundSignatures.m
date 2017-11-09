% GetConfoundSignatures
%-------------------------------------------------------------------------------
% Idea is to loop through a bunch of 'confounds' and determine their
% enrichment signatures -- can then be visualized alongside enrichment results
%-------------------------------------------------------------------------------

whatSpecies = 'human';

resultsTables = struct();

%-------------------------------------------------------------------------------
% Human -- decrease with distance
%-------------------------------------------------------------------------------
whatEdgeMeasure = 'distance';
onlyOnEdges = false;
correctDistance = false;
corrType = 'Spearman';
whatNull = 'randomGene';
numNulls = 100;

absType = 'neg';
[resultsTables.distance_neg,gScore] = EdgeEnrichment(whatEdgeMeasure,onlyOnEdges,...
                correctDistance,absType,corrType,whatNull,numNulls,whatSpecies)
absType = 'pos';
[resultsTables.distance_pos,gScore] = EdgeEnrichment(whatEdgeMeasure,onlyOnEdges,...
                correctDistance,absType,corrType,whatNull,numNulls,whatSpecies)

%===============================================================================
% Nodal correlations in cortex
%===============================================================================
structFilter = 'cortex';
corrType = 'Spearman';

% Highest variance:
enrichWhat = 'varExpression';
[resultsTables.varExpression,gScore] = NodeSimpleEnrichment(enrichWhat,structFilter,corrType,whatSpecies);

% Correlate with nodal degree:
enrichWhat = 'degree';
[resultsTables.degreeCorr,gScore] = NodeSimpleEnrichment(enrichWhat,structFilter,corrType,whatSpecies);

% Vary with dominant PC:
enrichWhat = 'genePC';
[resultsTables.PC1,gScore] = NodeSimpleEnrichment(enrichWhat,structFilter,corrType,whatSpecies);

%===============================================================================
% Now for some plotting:
%===============================================================================
theThreshold = 0.2;
PlotEnrichmentTables(resultsTables,theThreshold,whatSpecies);

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

PlotEnrichmentTables(resultsTablesTogether,theThreshold,whatSpecies);
