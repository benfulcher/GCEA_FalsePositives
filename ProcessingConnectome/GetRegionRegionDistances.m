function GetRegionRegionDistances()
% Retrieves region-region distances for all pairs of regions
% calculated in Suppl. Table 4 of the Oh et al. connectivity paper
% ------------------------------------------------------------------------------
% "Cartesian distances between the centers of mass of all the interconnected
% source and target region pairs for the 213 anatomical regions used in the
% linear model based Connectivity Matrix in figure 4a. Names of source regions
% (in rows) are shown in Column A, and names for target regions on both
% ipsilateral and contralateral hemispheres (in columns) are shown in Row 1. See
% Supplementary Table 1 for the corresponding full name and acronym of each
% region."
%
% This description is kind of confusing becasue there's no source/target
% distinction anymore, it's a symmetric matrix; rows and columns must be equivalent
% (ipsi and contra for both). I'm assuming it's not just for injections in a single
% hemisphere, but the all ipsi connections and all contra connections have been
% agglomerated...
%
% Also, have 295 regions/structures, whereas Fig. 4a has 213 regions/structures.
% Also, the caption for this supplement says "213 anatomical regions"... :/
%
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-05-29
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Loading...
% ------------------------------------------------------------------------------
% Read the Excel file:
InputFileName = 'nature13186-s5.xlsx';
fprintf(1,'Reading excel file %s\n',InputFileName);
[d_mat,LabelsRaw] = xlsread(InputFileName,'Distance');
fprintf(1,'Read\n');

% ------------------------------------------------------------------------------
% Organize the region labels from the excel file:
% ------------------------------------------------------------------------------
Labels = LabelsRaw(2:end,1); % get rows
NumLabels = length(Labels);

Reg_ipsi = regexp(Labels,'_ipsi');
IsIpsi = cellfun(@(x)~isempty(x),Reg_ipsi);
Reg_contra = regexp(Labels,'_contra');
IsContra = cellfun(@(x)~isempty(x),Reg_contra);

% Process each to be just the region name, and code each as either ipsi or contra
IpsiContra = zeros(NumLabels,1); % 1==ipsi, 2==contra
for i = 1:NumLabels
    if IsIpsi(i) & ~IsContra(i)
        Labels{i} = Labels{i}(1:Reg_ipsi{i}-1);
        IpsiContra(i) = 1;
    elseif ~IsIpsi(i) & IsContra(i)
        Labels{i} = Labels{i}(1:Reg_contra{i}-1);
        IpsiContra(i) = 2;
    else
        fprintf(1,'Weird one here...\n');
    end
end

% ------------------------------------------------------------------------------
% Match these regions to those used in the connectome
% ------------------------------------------------------------------------------
% Load Region acronyms from the Matlab file:
load('Mouse_Connectivity_Data.mat','regionAcronyms')

% This is a subset of regions represented in the distance matrix: need to match
UsedInConnectome = ismember(Labels,regionAcronyms);

% Filter based on regions used in the connectome:
Labels = Labels(UsedInConnectome);
IpsiContra = IpsiContra(UsedInConnectome);
d_mat = d_mat(UsedInConnectome,UsedInConnectome);

% Format into a 2x2 cell -- in each cell, the regions match the RegionStruct ordering
% {1,1: ipsi-ipsi}, {1,2: ipsi-contra}
% {2,1: contra-ipsi}, {2,2: contra-contra}
Dist_Matrix = cell(2,2);
for i = 1:2
    % Reordering of the labels for either ipsi (i==1) or contra (i==2)
    for j = 1:2
        % Do the reordering:
        MatchPoint_i = cellfun(@(x)find(strcmp(Labels(IpsiContra==i),x)),regionAcronyms);
        MatchPoint_j = cellfun(@(x)find(strcmp(Labels(IpsiContra==j),x)),regionAcronyms);

        Dist_Matrix{i,j} = d_mat(IpsiContra==i,IpsiContra==j);
        Dist_Matrix{i,j} = Dist_Matrix{i,j}(MatchPoint_i,MatchPoint_j);
    end
end

% So this Dist_Matrix is it -- rows/columns correspond to those in the RegionStruct

% ------------------------------------------------------------------------------
% Save back to file
% ------------------------------------------------------------------------------
% Dist_acronyms = Labels; % Acronym
% Dist_IpsiContra = IpsiContra; % Label for each entry: 1==ipsi, 2==contra
% Dist_Matrix = d_mat; % Distance matrix for all

save('Mouse_Connectivity_Data.mat','Dist_Matrix','-append')
fprintf(1,'Appended pairwise distance information to ''Mouse_Connectivity_Data.mat''\n');

end
