function [numSigMouse,numSigHuman] = LiteratureMatches(whatGOID,beVerbose)

if nargin < 1
    whatGOID = 8306;
end
if nargin < 2
    beVerbose = true;
end

%-------------------------------------------------------------------------------
% Tell them what it is:
params = GiveMeDefaultParams('mouse');
params.e.sizeFilter = [1,1e6];
GOTable = GiveMeGOData(params);
% Get the category of interest:
whatCategory = find(GOTable.GOID==whatGOID);
if isempty(whatCategory)
    error('No category found for ID %u using default mouse parameters',whatGOID);
end
fprintf(1,'GOGO %s\n',GOTable.GOName{whatCategory});

%-------------------------------------------------------------------------------
% Match to data:
load('LiteratureEnrichmentLoaded.mat','resultsTables','mouseOrHuman');

theAnalyses = fieldnames(resultsTables);
numAnalyses = length(theAnalyses);

numSigHuman = 0;
numSigMouse = 0;
for i = 1:numAnalyses
    if ismember(whatGOID,resultsTables.(theAnalyses{i}).GOID)
        if beVerbose
            fprintf(1,'Found in %s (%s)\n',theAnalyses{i},string(mouseOrHuman(i)));
        end
        theQuantity = resultsTables.(theAnalyses{i}).Properties.VariableNames{2};
        matchWhere = (resultsTables.(theAnalyses{i}).GOID==whatGOID);
        thePValue = resultsTables.(theAnalyses{i}).(theQuantity)(matchWhere);
        if beVerbose
            fprintf(1,'----%s = %g\n',theQuantity,thePValue);
        end
        if thePValue < 0.05
            if mouseOrHuman(i)=='mouse'
                numSigMouse = numSigMouse + 1;
            else
                numSigHuman = numSigHuman + 1;
            end
        end
    end
end
numSignificant = numSigMouse + numSigHuman;

end
