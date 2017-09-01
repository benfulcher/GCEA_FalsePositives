%-------------------------------------------------------------------------------
% Aim is to see evaluate the extent of overlap between enrichment results
% from a known nonspecific confound and some specific analysis
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Obtain the GOTables for degree:
GOTable_simple = NodeSimpleEnrichment('degree','all');
numSig = sum(GOTable_simple.pVal_corr < 0.05);
fprintf(1,'%u significant for degree (random gene nulls)\n',numSig);

GOTable_shuffle = NodeShuffleEnrichment('degree','all',200,'all');
numSig = sum(GOTable_shuffle.pValZ_corr < 0.05);
fprintf(1,'%u significant for degree (spatial nulls)\n',numSig);

% Obtain the GOTable for ranksum expression differences in isocortex:
GOTable_isocortex = NodeSimpleEnrichment('isocortex','all');
numSig = sum(GOTable_isocortex.pVal_corr < 0.05);
fprintf(1,'%u significant for cortex confound\n',numSig);

%-------------------------------------------------------------------------------
% Investigate overlaps between the different enrichment tables:
% Simplest is to see whether p-values are similar
