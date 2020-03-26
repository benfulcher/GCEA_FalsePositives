function TellMeLiteratureStory(whatGOID)
% Check inputs:
if nargin < 1
    whatGOID = 7158;
end

% Load info about this category:
params = GiveMeDefaultParams();
params.e.sizeFilter = [0,1e6];
GOTable = GiveMeGOData(params);
whatCategory = find(GOTable.GOID==whatGOID);
if isempty(whatCategory)
    error('No Category found with ID = %u',whatGOID);
end
fprintf(1,'GO Category: %s %s (%u)\n',GOTable.GOIDlabel{whatCategory},GOTable.GOName{whatCategory},GOTable.size(whatCategory));

%-------------------------------------------------------------------------------
% Find literature matches for this category:
load('LiteratureEnrichmentLoaded.mat','resultsTables','mouseOrHuman');
theStudyNames = fieldnames(resultsTables);
numEnrichmentStudies = length(theStudyNames);
for i = 1:numEnrichmentStudies
    % Look for the target in this study:
    if ismember(whatGOID,resultsTables.(theStudyNames{i}).GOID)
        fprintf(1,'Match for GOID %u in %s:     ',whatGOID,theStudyNames{i});
        here = (resultsTables.(theStudyNames{i}).GOID==whatGOID);
        fprintf(1,'pVal = %g\n',resultsTables.(theStudyNames{i}).pValCorr(here));
    end
end

end
