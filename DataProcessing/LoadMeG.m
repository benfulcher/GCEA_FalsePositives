function [GeneStruct,GeneExpData] = LoadMeG(removeDuplicates,doNormalize,energyOrDensity)
% Load in the gene data as G

%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 1 || isempty(removeDuplicates)
    removeDuplicates = true;
end
if nargin < 2 || isempty(doNormalize)
    doNormalize = 0;
end
if nargin < 3
    energyOrDensity = 'energy';
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Get the data:
%-------------------------------------------------------------------------------
allenDataPath = '/Users/benfulcher/GoogleDrive/Work/CompletedProjects/MouseConnectome/Code/';
fprintf(1,'Loading full gene data...');
load(fullfile(allenDataPath,'AllenGeneData_All.mat'),'GeneStruct','GeneExpData');
GeneExpData = GeneExpData.(energyOrDensity);
% 'RegionStruct',...
fprintf(1,' Done.\n');

%-------------------------------------------------------------------------------
% FURTHER PROCESSING:
%-------------------------------------------------------------------------------
if removeDuplicates
    [~,ia] = unique([GeneStruct.gene_entrez_id]);
    fprintf(1,'Filtering down genes from %u to a unique set of %u\n',...
                        length([GeneStruct.gene_entrez_id]),length(ia));
    GeneExpData = GeneExpData(:,ia);
    GeneStruct = GeneStruct(ia);
end
%-------------------------------------------------------------------------------
if doNormalize
    GeneExpData = BF_NormalizeMatrix(GeneExpData,'mixedSigmoid');
    fprintf(1,'Normalized using mixed sigmoid\n');
end

end
