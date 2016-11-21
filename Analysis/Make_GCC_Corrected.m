function MakeGCC()
%-------------------------------------------------------------------------------
% Quick function to try to extract a smaller file containing GCC scores
% for all genes, only at connections
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Set Parameters:
whatGeneData = 'raw'; % 'norm', 'raw', 'zscore'
energyOrDensity = 'energy'; % 'energy', 'density'

%-------------------------------------------------------------------------------
% Load in gene data

%-------------------------------------------------------------------------------
% Get the elements where links of binary A at 0.05 exist:
C = load('Mouse_Connectivity_Data.mat'); % C stands for connectome data
% Get the adjacency matrix and process it:
A = logical(GiveMeAdj(C,'binary','ipsi',0,0.05));

if removeDeadEnd
    k_out = sum(A,2); % out degree
    deadEndNodes = (k_out==0);
    A = A(~deadEndNodes,~deadEndNodes);
    ggBlockscorr = ggBlockscorr(:,~deadEndNodes,~deadEndNodes);
    fprintf(1,'Removed %u dead end nodes\n',sum(deadEndNodes));
end

%-------------------------------------------------------------------------------
numLinks = sum(A(:));
ggCorr_A_p005 = zeros(numGenes,numLinks);
for i = 1:numGenes
    GCC_i = squeeze(ggBlockscorr(i,:,:));
    ggCorr_A_p005(i,:) = GCC_i(A);
end


if removeDeadEnd
    fileName = 'GCC_A_p005_deadEnd.mat';
else
    fileName = 'GCC_A_p005.mat';
end

save(fileName,'ggCorr_A_p005','theGeneStruct','whatGeneData','energyOrDensity',...
                'theCorrection','theFitType','baselineLinks','removeDeadEnd');

end
