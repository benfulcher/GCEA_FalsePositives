% Add paths required for the project (ignoring hidden, including version control)
files = dir;
directories = files([files.isdir]);
directories(strmatch('.',{files([files.isdir]).name})) = []; % remove hidden
paths = arrayfun(@(x)fullfile(directories(x).folder,directories(x).name),1:length(directories),'UniformOutput',false);
for j = 1:length(paths)
    addpath(genpath(paths{j}))
end

%===============================================================================
% Add path references to dependencies:
%===============================================================================
% Matlab gene enrichment toolbox (github):
addpath('/Users/benfulcher/GoogleDrive/Work/CodeToolboxes/MatlabEnrichment/')
addpath('/Users/benfulcher/GoogleDrive/Work/CodeToolboxes/MatlabmySQL/')
