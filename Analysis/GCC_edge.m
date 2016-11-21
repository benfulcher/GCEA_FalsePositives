% function GCC_edge()

%-------------------------------------------------------------------------------
% Set Parameters:
pThreshold = 0.05; % for connectivity data

%-------------------------------------------------------------------------------
% Load basic connectivity data:
C = load('Mouse_Connectivity_Data.mat'); % C stands for connectome data

%-------------------------------------------------------------------------------
% Get the adjacency matrix and process it:
A_weight = GiveMeAdj(C,'zero','ipsi',0,pThreshold);
A_bin = (A_weight > 0);
% Filter out dead-end nodes:
k_out = sum(A_bin,2); % out degree
deadEndNodes = (k_out==0);
A_bin = A_bin(~deadEndNodes,~deadEndNodes);
A_weight = A_weight(~deadEndNodes,~deadEndNodes);

%-------------------------------------------------------------------------------
% Load the GCC scores from file
load('GCC_A_p005_deadEnd.mat','ggCorr_A_p005','theGeneStruct');
numGenes = size(ggCorr_A_p005,1);

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Compute all edge measures
edgeMeasures = GiveMeAllEdgeMeasures(A_bin,A_weight);
numStats = length(fieldnames(edgeMeasures));

%-------------------------------------------------------------------------------
% Ok, so now we can run stats on the corrected ggBlockscorr (GCC scores)
gScore = zeros(numGenes,numStats);
scoreNames = fieldnames(edgeMeasures);
q = zeros(sum(A_bin(:)),numStats);
for i = 1:numStats
    if ismember(scoreNames{i},{'edgeBet','G','signalTraffic'})
        fprintf(1,'Doing a log10 transform for %s\n',scoreNames{i});
        q(:,i) = log10(edgeMeasures.(scoreNames{i})(A_bin));
    else
        q(:,i) = edgeMeasures.(scoreNames{i})(A_bin);
    end
end
for i = 1:numGenes
    GCC_i = ggCorr_A_p005(i,:)';
    % Keep only good values
    goodData = ~isnan(GCC_i);
    GCC_i = GCC_i(goodData);
    for j = 1:numStats
        gScore(i,j) = corr(q(goodData,j),GCC_i); % log10 edge betweenness
    end
end

%-------------------------------------------------------------------------------
% Save quickly to .mat file:
entrezIDs = [theGeneStruct.filtered.gene_entrez_id];
save('gScore.mat','gScore','entrezIDs','scoreNames')

%-------------------------------------------------------------------------------
% Save each result to a separate ErmineJ file:
for i = 1:length(scoreNames)
    fileName = sprintf('%s',scoreNames{i});
    writeErmineJFile(fileName,gScore(:,i),entrezIDs,scoreNames{i});
end

% end
