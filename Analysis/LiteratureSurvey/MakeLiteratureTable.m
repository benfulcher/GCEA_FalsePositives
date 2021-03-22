function LiteratureGOTable = MakeLiteratureTable(whatSpecies,pValCorrThreshold)
%-------------------------------------------------------------------------------
% Information about the most implicated GO categories in the literature
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Filter by a particular species:
if nargin < 1
    whatSpecies = 'human';
    warning('human by default')
end
if nargin < 2
    pValCorrThreshold = 0.05;
end
%-------------------------------------------------------------------------------

% Load literature annotations:
load('LiteratureEnrichmentLoaded.mat','resultsTables','mouseOrHuman');

% Give some basic info

% Filter on species:
isSpeciesOfinterest = (mouseOrHuman==whatSpecies);
fNames = fieldnames(resultsTables);
resultsTables = rmfield(resultsTables,fNames(~isSpeciesOfinterest));

% Thresholds so that (i) only annotations with below-threshold significance are included, and
% (ii) only studies containing such below-threshold categories are filtered
% -> resultsTablesTh
resultsTablesTh = resultsTables;
resultsTablesTh = structfun(@(x)x(x.pValCorr<=pValCorrThreshold,:),resultsTablesTh,'UniformOutput',false);
emptyIndex = find(structfun(@(x)(height(x)==0),resultsTablesTh));
numEmpty = length(emptyIndex);
allFields = fieldnames(resultsTablesTh);
emptyFields = allFields(emptyIndex);
for i = 1:numEmpty
    resultsTablesTh = rmfield(resultsTablesTh,emptyFields{i});
end

% Convert to a table of GO categories, annotated by a list of studies
allFlaggedGOIDs = structfun(@(x)x.GOID,resultsTablesTh,'UniformOutput',false);
allFlaggedGOIDs = struct2cell(allFlaggedGOIDs);
allFlaggedGOIDs = unique(vertcat(allFlaggedGOIDs{:}));

%-------------------------------------------------------------------------------
% Now form a table:
numGOIDs = length(allFlaggedGOIDs);
studyList = cell(numGOIDs,1);
allFields = fieldnames(resultsTablesTh);
for i = 1:numGOIDs
    isListed = structfun(@(x)ismember(allFlaggedGOIDs(i),x.GOID),resultsTablesTh);
    studyList{i} = allFields(isListed);
end
GOID = allFlaggedGOIDs;
numStudies = cellfun(@length,studyList);
studyListMessyString = cellfun(@(x)BF_cat(x),studyList,'UniformOutput',false);
LiteratureGOTable = table(GOID,numStudies,studyListMessyString,studyList);
[~,ix] = sort(numStudies,'descend');
LiteratureGOTable = LiteratureGOTable(ix,:);

end
