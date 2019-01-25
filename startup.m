% Add paths required for the project (ignoring hidden, including version control)
files = dir;
directories = files([files.isdir]);
directories(strmatch('.',{files([files.isdir]).name})) = []; % remove hidden
paths = arrayfun(@(x)fullfile(directories(x).folder,directories(x).name),1:length(directories),'UniformOutput',false);
for j = 1:length(paths)
    addpath(genpath(paths{j}))
end

%===============================================================================
% Add path to dependencies:

fprintf(1,'Adding dependencies for external toolboxes:\n');

fprintf(1,'GeneEnrichment for Matlab\n');
pathToToolbox = GiveMeFile('EnrichmentToolbox');
addpath(pathToToolbox)
hereNow = pwd;
cd(pathToToolbox)
startup;
cd(hereNow);

fprintf(1,'mySQL for Matlab\n');
try
    addpath('/Users/benfulcher/DropboxSydneyUni/CodeToolboxes/MatlabmySQL/')
catch
    fprintf(1,'Failed\n');
end
