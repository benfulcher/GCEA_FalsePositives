function [W_rect,sourceRegions,targetRegions] = ImportCorticalConnectivityWeights()
% In the mouse dataset, 31 of the 43 cortical regions were injected with at least
% half the injection volume in at least one experiment.
% We obtained point estimates and credibility intervals for 2,627 out of
% 31 × 86 − 31 = 2,635 (99.7%) possible pair-wise con- nections;
% see S3 Table for full connectivity matrix.
% Self connections were not measured, the remaining 8 connections (i,j) could
% not be measured as the target region j received part of the injection for all
% experiments where i was injected
% (from Ypma and Bullmore, 2016)

%-------------------------------------------------------------------------------
%% First the connectivity weights:
[~, ~, raw0] = xlsread('Ypma_cortical_connectome_journal.pcbi.1005104.s003.xls',...
                                'log NCD estimates');
raw = raw0(3:end,45:end);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
W_rect = reshape([raw{:}],size(raw));
% Convert from log10 to real numbers:
W_rect = arrayfun(@(x)10^x,W_rect);

%-------------------------------------------------------------------------------
%% Next the region names:
sourceRegions = raw0(3:end,1); % has to be this way around because left hemisphere are columns
targetRegions = raw0(2,45:end)';

end
