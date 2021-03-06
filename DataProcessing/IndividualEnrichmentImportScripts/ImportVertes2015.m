function FilteredTable = ImportVertes2015(whichSheet)
if nargin < 1
    whichSheet = 'PLS1 pos';
end
%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/benfulcher/GoogleDrive/Work/CurrentProjects/GeneExpressionEnrichment/DataSets/Vertes-rstb20150362supp1.xlsx
%    Worksheet: PLS1 pos
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Auto-generated by MATLAB on 2017/10/27 16:54:38

%% Import the data
[~, ~, raw] = xlsread('Vertes-rstb20150362supp1.xlsx',whichSheet);
raw = raw(3:end,3:8);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
stringVectors = string(raw(:,[1,2,5]));
stringVectors(ismissing(stringVectors)) = '';
raw = raw(:,[3,4,6]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Create table
ResultsTable = table;

%% Allocate imported array to column variable names
ResultsTable.Description = stringVectors(:,1);
ResultsTable.GOterm = stringVectors(:,2);
ResultsTable.FDRqvalue = data(:,1);
ResultsTable.Pvalue = data(:,2);
ResultsTable.EnrichmentNBnb = stringVectors(:,3);
ResultsTable.exclude = data(:,3);

%-------------------------------------------------------------------------------
% Now just take the necessary columns
isEmpty = cellfun(@isempty,ResultsTable.GOterm);
ResultsTable = ResultsTable(~isEmpty,:);
pValCorr = ResultsTable.FDRqvalue;
GOtoNumber = @(x)str2num(x(4:10));
GOID = cellfun(GOtoNumber,ResultsTable.GOterm);
FilteredTable = sortrows(table(GOID,pValCorr),'pValCorr');

end
