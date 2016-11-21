function G = LoadMeG()
% load in the gene data as G

allenDataPath = '/Users/benfulcher/GoogleDrive/Work/CompletedProjects/MouseConnectome/Code/';
fprintf(1,'Loading full gene data...');
G = load(fullfile(allenDataPath,'AllenGeneData_All.mat'),'RegionStruct','GeneStruct','GeneExpData');
fprintf(1,' Done.\n');

end
