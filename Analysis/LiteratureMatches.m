function LiteratureMatches(whatGOID)

if nargin < 1
    whatGOID = 8306;
end

%-------------------------------------------------------------------------------
% Tell them what it is:
params = GiveMeDefaultParams('mouse');
GOTable = GiveMeGOData(params);
% Get the category of interest:
whatCategory = find(GOTable.GOID==whatGOID);
fprintf(1,'GOGO %s\n',GOTable.GOName{whatCategory});

%-------------------------------------------------------------------------------
% Match to data:
load('LiteratureEnrichmentLoaded.mat','resultsTables','mouseOrHuman');

theAnalyses = fieldnames(resultsTables);
numAnalyses = length(theAnalyses);

for i = 1:numAnalyses
    if ismember(whatGOID,resultsTables.(theAnalyses{i}).GOID)
        fprintf(1,'Found in %s\n',theAnalyses{i});
        theQuantity = resultsTables.(theAnalyses{i}).Properties.VariableNames{2};
        matchWhere = (resultsTables.(theAnalyses{i}).GOID==whatGOID);
        fprintf(1,'----%s = %g\n',theQuantity,...
                    resultsTables.(theAnalyses{i}).(theQuantity)(matchWhere));
    end
end

end
