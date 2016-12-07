function [GeneStruct,GeneExpData] = LoadMeG_old(removeDuplicates,normalizationSettings,energyOrDensity)
% Load in the gene data as G (FROM OLD RETRIEVAL)

%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 1 || isempty(removeDuplicates)
    removeDuplicates = true;
end
if nargin < 2 || isempty(normalizationSettings)
    normalizationSettings = {};
    % In the format: {howToNormalizeGenesAcrossRegions,howToNormalizeRegionsAcrossGenes}
end
if nargin < 3
    energyOrDensity = 'energy';
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Get the OLD data:
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
%-------------------------------------------------------------------------------
if removeDuplicates
    [~,ia] = unique([GeneStruct.gene_entrez_id]);
    fprintf(1,'Filtering down genes from %u to a unique set of %u\n',...
                        length([GeneStruct.gene_entrez_id]),length(ia));
    GeneExpData = GeneExpData(:,ia);
    GeneStruct = GeneStruct(ia);
end
if ~isempty(normalizationSettings)
    GeneExpData = BF_NormalizeMatrix(GeneExpData,normalizationSettings{1});
    fprintf(1,'1. Normalized expression for each gene using %s\n',normalizationSettings{1});
    GeneExpData = BF_NormalizeMatrix(GeneExpData',normalizationSettings{2})';
    fprintf(1,'2. Normalized expression across each brain region using %s\n',normalizationSettings{2});
end

end
