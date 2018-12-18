function resultsTablesSpecies = LiteratureLook(whatSpecies,theThreshold,doPlot)
% Unifies across literature enrichment results and displays in a graphical table
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    whatSpecies = 'human'; % 'mouse','human'
end
if nargin < 2
    % Threshold for displaying a category as "significant"
    theThreshold = 0.2;
end
if nargin < 3
    doPlot = true;
end
%-------------------------------------------------------------------------------

% Load in data on enrichment results across studies:
% (saved from ImportLiteratureEnrichment)
load('LiteratureEnrichmentLoaded.mat','resultsTables','mouseOrHuman');

% Filter just to tables involving the specified species:
allTableNames = fieldnames(resultsTables);
theSpeciesTables = find(mouseOrHuman==whatSpecies);
speciesTableNames = allTableNames(theSpeciesTables);
numSpecies = length(theSpeciesTables);
resultsTablesSpecies = struct;
for i = 1:numSpecies
    theName = speciesTableNames{i};
    resultsTablesSpecies.(theName) = resultsTables.(theName);
end
fprintf(1,'%u GO enrichment datasets involving %s\n',numSpecies,whatSpecies);

if doPlot
    PlotEnrichmentTables(resultsTablesSpecies,theThreshold);
    title(whatSpecies);
end

end
