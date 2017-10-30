function FilteredResults = ImportWhitaker(whatSheet)
%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/benfulcher/GoogleDrive/Work/CurrentProjects/GeneExpressionEnrichment/DataSets/Whitaker/WhitakerReformatted.xlsx
%    Worksheet: Complete_PLS2pos
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Auto-generated by MATLAB on 2017/10/28 11:49:52

if nargin < 1
    whatSheet = 'Complete_PLS2pos';
end

%% Import the data
[~, ~, raw] = xlsread('/Users/benfulcher/GoogleDrive/Work/CurrentProjects/GeneExpressionEnrichment/DataSets/Whitaker/WhitakerReformatted.xlsx',whatSheet);
raw = raw(2:end,:);
stringVectors = string(raw(:,[1,2,5]));
stringVectors(ismissing(stringVectors)) = '';
raw = raw(:,[3,4,6]);

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Create table
ResultsTable = table;

%% Allocate imported array to column variable names
ResultsTable.GOtermID = stringVectors(:,1);
ResultsTable.GOterm = stringVectors(:,2);
ResultsTable.rawPvalue = data(:,1);
ResultsTable.FDRqvalue = data(:,2);
ResultsTable.EnrichmentNBnb = stringVectors(:,3);
ResultsTable.isTooBig = data(:,3);

%-------------------------------------------------------------------------------
% Now just take the necessary columns
isEmpty = cellfun(@isempty,ResultsTable.GOtermID);
ResultsTable = ResultsTable(~isEmpty,:);
pValCorr = ResultsTable.FDRqvalue;
GOtoNumber = @(x)str2num(x(4:end));
GOID = cellfun(GOtoNumber,ResultsTable.GOtermID);
FilteredResults = sortrows(table(GOID,pValCorr),'pValCorr');

end