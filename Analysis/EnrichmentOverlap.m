%-------------------------------------------------------------------------------
% Aim is to see evaluate the extent of overlap between enrichment results
% from a known nonspecific confound and some specific analysis
%-------------------------------------------------------------------------------
% If we produce confound lists, we could have a principled (semi-automatic?) way
% to quantify whether a specific analysis may be driven by a given confound by
% the similarity in category scores...
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Obtain the GOTables for degree:
GOTable_simple = NodeSimpleEnrichment('degree','all');
numSig = sum(GOTable_simple.pValCorr < 0.05);
fprintf(1,'%u significant for degree (random gene nulls)\n',numSig);

GOTable_shuffle = NodeShuffleEnrichment('degree','all',200,'all');
numSig = sum(GOTable_shuffle.pValZCorr < 0.05);
fprintf(1,'%u significant for degree (spatial nulls)\n',numSig);

% Obtain the GOTable for ranksum expression differences in isocortex:
GOTable_isocortex = NodeSimpleEnrichment('isocortex','all');
numSig = sum(GOTable_isocortex.pValCorr < 0.05);
fprintf(1,'%u significant for cortex confound\n',numSig);

%-------------------------------------------------------------------------------
% Investigate overlaps between the different enrichment tables:
% Simplest is to see whether p-values are similar
PlotGOScoreScatter(GOTable_simple,GOTable_isocortex);
xlabel('GO category score (mean Spearman correlation with degree)')
ylabel('GO category score (-log10 p-value ranksum test isocortex)')
