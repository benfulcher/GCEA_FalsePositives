function LoadConnectomeData()
% ------------------------------------------------------------------------------
% ImportData
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-05-09 -- imports the Connectome data from Excel
% ------------------------------------------------------------------------------

InputFileName = 'Suppl_Table_3_nature13186-s4.xlsx';

SheetNames = {'W_ipsi','PValue_ipsi','W_contra','PValue_contra'};

% ------------------------------------------------------------------------------
% Import Excel data from all four sheets:

for WhatSheet = 1:length(SheetNames)
    [Adj{WhatSheet},RegionAcronyms{WhatSheet}] = xlsread(InputFileName,SheetNames{WhatSheet});

    fprintf(1,'Data imported from Excel file ''%s:%s''!\n', ...
                        InputFileName,SheetNames{WhatSheet});

    RowLabels = RegionAcronyms{WhatSheet}(2:end,1);
    ColumnLabels = RegionAcronyms{WhatSheet}(1,2:end)';
    NumRegions(WhatSheet) = size(Adj{WhatSheet},1);

    Match = arrayfun(@(x)strcmp(RowLabels(x),ColumnLabels(x)),1:NumRegions);

    if any(Match==0)
        error('Column and row labels do not match... :-/');
    else
        fprintf(1,'All row and column names match\n');
    end
    RegionAcronyms{WhatSheet} = RowLabels;
    clear('RowLabels','ColumnLabels','Match')

end
clear('WhatSheet')


% ------------------------------------------------------------------------------
% Assimilate data
% ------------------------------------------------------------------------------
if any(NumRegions~=NumRegions(1))
    error('Number of regions don''t match across Excel sheets');
else
    NumRegions = NumRegions(1);
end

MatchingLabels = zeros(length(SheetNames),length(SheetNames));
for i = 1:length(SheetNames)-1
    for j = i:length(SheetNames)
        MatchingLabels(i,j) = any(arrayfun(@(x)strcmp(RegionAcronyms{i}(x), ...
                                RegionAcronyms{j}(x)),1:length(RegionAcronyms{1}))~=1);
    end
end

if any(MatchingLabels(:)~=0)
    error('Labels contra/ipsi don''t match');
else
    fprintf(1,'Labels match in all four sheets\n');
    RegionAcronyms = RegionAcronyms{1};
end
clear('MatchingLabels')


% ------------------------------------------------------------------------------
% Reorder according to Oh default major region ordering
% ------------------------------------------------------------------------------
% Get the RegionStruct to match to:
fprintf(1,'Loading RegionStruct to match regions to...\n');
theMatFile = 'AllenGeneData_Coronal.mat';
fprintf(1,'************* Loading Oh-ordered Region structure information from %s\n',theMatFile);
load(theMatFile,'RegionStruct');
% Reorder Adj, and RegionAcronyms:
perm_c_reg = cellfun(@(x)find(strcmp(x,RegionAcronyms)),{RegionStruct.acronym});
RegionAcronyms = RegionAcronyms(perm_c_reg);
for i = 1:length(Adj)
    Adj{i} = Adj{i}(perm_c_reg,perm_c_reg);
end
fprintf(1,'Reordered to match RegionStruct!\n');

% ------------------------------------------------------------------------------
% Normalize
% ------------------------------------------------------------------------------
fprintf(1,'Normalizing data....\n');
fprintf(1,'Using a iqr-sigmoid to normalize non-zero and non-NaN entries as if they''re a vector...\n');
Adj_norm = cell(2,1);
for i = 1:4
    r = (~isnan(Adj{i}) & Adj{i} > 0);
    allWeights = Adj{i}(r); % positive, non-NaN weights
    allWeights_norm = BF_NormalizeMatrix(allWeights,'MixedSigmoid');
    Adj_norm{i} = zeros(size(Adj{i}));
    Adj_norm{i}(r) = allWeights_norm;
end
clear r allWeights allWeights_norm
fprintf(1,'Done\n');

% ------------------------------------------------------------------------------
% Cluster
% ------------------------------------------------------------------------------
fprintf(1,'Clustering by Euclidean distances...');
% Cluster rows by Euclidean distances
Cluster_ix = cell(4,1);
Cluster_norm_ix = cell(4,1);
for i = 1:4
    % Cluster raw values (Cluster_ix):
    ix = Cluster_Reorder(Adj{i},'euclidean','average');
    Cluster_ix{i} = ix; % Save the permutation
    Adj_cl{i} = Adj{i}(ix,ix);

    % Cluster normalized values (Cluster_norm_ix):
    ix = Cluster_Reorder(Adj_norm{i},'euclidean','average');
    Cluster_norm_ix{i} = ix; % Save the permutation
    Adj_norm_cl{i} = Adj_norm{i}(ix,ix);
end
clear ix
fprintf(1,'Done\n');

% ------------------------------------------------------------------------------
% Separate into p-value and adjacency matrix components
% ------------------------------------------------------------------------------
Conn_W = Adj([1,3]);
Conn_p = Adj([2,4]);
Conn_W_norm = Adj_norm([1,3]);
Conn_p_norm = Adj_norm([2,4]);
Conn_W_norm_cl = Adj_norm_cl([1,3]);
Conn_p_norm_cl = Adj_norm_cl([2,4]);
Conn_W_cl = Adj_cl([1,3]);
Conn_p_cl = Adj_cl([2,4]);

% Let's keep the Adj stuffs for old routines that use it...:
% clear Adj Adj_norm Adj_cl Adj_norm_cl

% Clustering permutations:
Cluster_ix_W = Cluster_ix([1,3]);
Cluster_ix_p = Cluster_ix([2,4]);
Cluster_norm_ix_W = Cluster_norm_ix([1,3]);
Cluster_norm_ix_p = Cluster_norm_ix([2,4]);

clear Cluster_ix Cluster_norm_ix

% ------------------------------------------------------------------------------
% Threshold
% ------------------------------------------------------------------------------
% Threshold the connections at a given p-value; ignore others

p_th = 0.05;
Conn_W_th = cell(2,1); % ipsi,contra
Conn_W_norm_th = cell(2,1); % ipsi,contra

for i = 1:2 % ipsi,contra
    % Raw connection weights, thresholded on p-values
    Conn_W_th{i} = ConnectivityAtPThreshold(Conn_W{i},Conn_p{i},p_th,NaN);

    % Normalized connection weights, thresholded on p-values
    % Threshold first, because normalization transformation will change depending
    % on the distribution of remaining weights:
    r = (~isnan(Conn_W_th{i}) & Conn_W_th{i} > 0);
    allWeights = Conn_W_th{i}(r); % positive, non-NaN weights
    allWeights_norm = BF_NormalizeMatrix(allWeights,'MixedSigmoid');
    Conn_W_th_norm{i} = zeros(size(Conn_W_th{i}));
    Conn_W_th_norm{i}(r) = allWeights_norm;

    % Threshold pre-normalized threshold weights:
    Conn_W_norm_th{i} = ConnectivityAtPThreshold(Conn_W_norm{i},Conn_p{i},p_th,NaN);
end
clear r allWeights allWeights_norm

% ------------------------------------------------------------------------------
% Make p-value thresholded, zeroed versions
% ------------------------------------------------------------------------------
Conn_W_zero = Conn_W_th;
Conn_W_norm_zero = Conn_W_th_norm;
for i = 1:2
    Conn_W_zero{i}(isnan(Conn_W_zero{i})) = 0;
    Conn_W_norm_zero{i}(isnan(Conn_W_norm_zero{i})) = 0;
end

% ------------------------------------------------------------------------------
% Make binary versions
% ------------------------------------------------------------------------------
Conn_bin = arrayfun(@(x) Conn_W_zero{x} > 0,1:2,'UniformOutput',0);

% ------------------------------------------------------------------------------
% Make p-value weights
% Ben Fulcher, 2014-08-22
% ------------------------------------------------------------------------------
% Defines weight as 1-p (where missing links have weight 0).
makeComplement = @(x) 1-x;
Conn_pweight = arrayfun(@(x) makeComplement(makeNaN1(Conn_p{x})),1:2,'UniformOutput',0);
clear makeComplement

% ------------------------------------------------------------------------------
% Make combination p-value weighted weights!
% Ben Fulcher, 2014-08-22
% ------------------------------------------------------------------------------
% high weight is product of p-value with the weight value.
% Probably not theoretically sound, but it's workable...
Conn_W_pweighted = arrayfun(@(x) Conn_pweight{x}.*Conn_W{x},1:2,'UniformOutput',0);

% ------------------------------------------------------------------------------
% Save to file
% ------------------------------------------------------------------------------
clear i j
save('Mouse_Connectivity_Data.mat')
fprintf(1,'Data saved to %s!\n','Mouse_Connectivity_Data.mat');

% ------------------------------------------------------------------------------
% Add major region labels
% ------------------------------------------------------------------------------
LabelMajorBrainRegions

% ------------------------------------------------------------------------------
% Add region-region Euclidean distances
% ------------------------------------------------------------------------------
GetRegionRegionDistances

%-------------------------------------------------------------------------------
function x = makeNaN1(x)
    x(isnan(x)) = 1;
end

end
