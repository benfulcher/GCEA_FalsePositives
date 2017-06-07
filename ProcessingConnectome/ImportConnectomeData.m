function ImportConnectomeData()
% ------------------------------------------------------------------------------
% ImportData
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-05-09 -- imports the Connectome data from Excel
% ------------------------------------------------------------------------------

InputFileName = 'Suppl_Table_3_nature13186-s4.xlsx';

SheetNames = {'W_ipsi','PValue_ipsi','W_contra','PValue_contra'};

% ------------------------------------------------------------------------------
% Import Excel data from all four sheets:

for whatSheet = 1:length(SheetNames)
    [Adj{whatSheet},regionAcronyms{whatSheet}] = xlsread(InputFileName,SheetNames{whatSheet});

    fprintf(1,'Data imported from Excel file ''%s:%s''!\n', ...
                        InputFileName,SheetNames{whatSheet});

    rowLabels = regionAcronyms{whatSheet}(2:end,1);
    columnLabels = regionAcronyms{whatSheet}(1,2:end)';
    numRegions(whatSheet) = size(Adj{whatSheet},1);

    match = arrayfun(@(x)strcmp(rowLabels(x),columnLabels(x)),1:numRegions);

    if any(match==0)
        error('Column and row labels do not match... :-/');
    else
        fprintf(1,'All row and column names match\n');
    end
    regionAcronyms{whatSheet} = rowLabels;
end

% ------------------------------------------------------------------------------
% Assimilate data
% ------------------------------------------------------------------------------
if any(numRegions~=numRegions(1))
    error('Number of regions don''t match across Excel sheets');
else
    numRegions = numRegions(1);
end

matchingLabels = zeros(length(SheetNames),length(SheetNames));
for i = 1:length(SheetNames)-1
    for j = i:length(SheetNames)
        matchingLabels(i,j) = any(arrayfun(@(x)strcmp(regionAcronyms{i}(x), ...
                                regionAcronyms{j}(x)),1:length(regionAcronyms{1}))~=1);
    end
end

if any(matchingLabels(:)~=0)
    error('Labels contra/ipsi don''t match');
else
    fprintf(1,'Labels match in all four sheets\n');
    regionAcronyms = regionAcronyms{1};
end
clear('matchingLabels')


% ------------------------------------------------------------------------------
% Reorder according to region ordering in gene expression data
% ------------------------------------------------------------------------------
% Get the RegionStruct to match to:
dataFile = '/Users/benfulcher/GoogleDrive/Work/CurrentProjects/CellTypesMouse/Code/AllenGeneDataset_19419.mat';
load(dataFile,'structInfo');

% Reorder Adj, and regionAcronyms:
[~,iG,ix] = intersect(structInfo.acronym,regionAcronyms,'stable');
if length(iG)~=length(ix)
    error('Error matching structures');
end
regionAcronyms = regionAcronyms(ix);
for i = 1:length(Adj)
    Adj{i} = Adj{i}(ix,ix);
end
fprintf(1,'Reordered to match gene data!\n');

% ------------------------------------------------------------------------------
% Separate into p-value and adjacency matrix components
% ------------------------------------------------------------------------------
Conn_W = Adj([1,3]);
Conn_p = Adj([2,4]);

%-------------------------------------------------------------------------------
% Save back to file:
%-------------------------------------------------------------------------------
save(fullfile('DataOutputs','Mouse_Connectivity_Data.mat'),'Conn_W','Conn_p','regionAcronyms')

% ------------------------------------------------------------------------------
% Add region-region Euclidean distances
% ------------------------------------------------------------------------------
GetRegionRegionDistances

%-------------------------------------------------------------------------------
function x = makeNaN1(x)
    x(isnan(x)) = 1;
end

end
