function G = LoadMeG(removeDuplicates)
% Load in the gene data as G
if nargin < 1 || isempty(removeDuplicates)
    removeDuplicates = true;
end

allenDataPath = '/Users/benfulcher/GoogleDrive/Work/CompletedProjects/MouseConnectome/Code/';
fprintf(1,'Loading full gene data...');
G = load(fullfile(allenDataPath,'AllenGeneData_All.mat'),'RegionStruct',...
                                            'GeneStruct','GeneExpData');
fprintf(1,' Done.\n');

if removeDuplicates
    [~,ia] = unique([G.GeneStruct.gene_entrez_id]);
    fprintf(1,'Filtering down genes from %u to a unique set of %u\n',...
                        length([G.GeneStruct.gene_entrez_id]),length(ia));
    G.GeneExpData.energy = G.GeneExpData.energy(:,ia);
    G.GeneExpData.density = G.GeneExpData.density(:,ia);
    G.GeneStruct = G.GeneStruct(ia);
end

end
