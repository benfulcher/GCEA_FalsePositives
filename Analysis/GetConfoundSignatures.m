% GetConfoundSignatures
%-------------------------------------------------------------------------------
% Idea is to loop through a bunch of 'confounds' and determine their
% enrichment signatures
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Human -- decrease with distance
%-------------------------------------------------------------------------------
whatEdgeMeasure = 'distance';
onlyOnEdges = false;
correctDistance = false;
absType = 'neg';
corrType = 'Spearman';
whatNull = 'randomGene';
numNulls = 100;
whatSpecies = 'human';
[GOTable_distance,gScore] = EdgeEnrichment(whatEdgeMeasure,onlyOnEdges,...
                correctDistance,absType,corrType,whatNull,numNulls,whatSpecies)


%-------------------------------------------------------------------------------
% Human -- increase with degree
%-------------------------------------------------------------------------------
enrichWhat = 'degree';
structFilter = 'cortex';
whatCorr = 'Spearman';
whatSpecies = 'human';
[GOTable_k,gScore] = NodeSimpleEnrichment(enrichWhat,structFilter,whatCorr,whatSpecies);

%-------------------------------------------------------------------------------
% Human -- vary with PC
%-------------------------------------------------------------------------------
enrichWhat = 'genePC';
structFilter = 'cortex';
whatSpecies = 'human';
whatCorr = 'Spearman';
[GOTable_k,gScore] = NodeSimpleEnrichment(enrichWhat,structFilter,whatCorr,whatSpecies);
