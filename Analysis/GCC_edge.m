% function GCC_edge()

%-------------------------------------------------------------------------------
% Set parameters:
pThreshold = 0.05; % for connectivity data
whatCorr = 'Spearman'; % what correlation metric

%-------------------------------------------------------------------------------
% Load connectivity data and process adjacency matrix
C = load('Mouse_Connectivity_Data.mat'); % C stands for connectome data
A_weight = GiveMeAdj(C,'zero','ipsi',0,pThreshold);
A_bin = (A_weight > 0);
% Filter out dead-end nodes:
k_out = sum(A_bin,2); % out degree
deadEndNodes = (k_out==0);
A_bin = A_bin(~deadEndNodes,~deadEndNodes);
A_weight = A_weight(~deadEndNodes,~deadEndNodes);

%-------------------------------------------------------------------------------
% Get gene expression data:
G = LoadMeG(true);
GData = G.GeneExpData.(energyOrDensity);
numGenes = size(GData,2);

%-------------------------------------------------------------------------------
% Compute all edge measures
edgeMeasures = GiveMeAllEdgeMeasures(A_bin,A_weight);
edgeMeasureNames = fieldnames(edgeMeasures);
numStats = length(edgeMeasureNames);

% Convert edge measures to matrix with a row for each connection:
q = zeros(sum(A_bin(:)),numStats);
for i = 1:numStats
    q(:,i) = edgeMeasures.(edgeMeasureNames{i})(A_bin);
    % Transforms necessary for Pearson, e.g.,
    % if ismember(edgeMeasureNames{i},{'edgeBet','G','signalTraffic'})
    %     fprintf(1,'Doing a log10 transform for %s\n',edgeMeasureNames{i});
    %     q(:,i) = log10(edgeMeasures.(edgeMeasureNames{i})(A_bin));
    % else
    %     q(:,i) = edgeMeasures.(edgeMeasureNames{i})(A_bin);
    % end
end

%-------------------------------------------------------------------------------
% Ok, so now we can find correlations to GCC scores across genes
gScore = zeros(numGenes,numStats);
fprintf(1,'Looping over %u genes, %u edge-based statistics\n',numGenes,numStats);
for i = 1:numGenes
    g = GData(:,i);
    GCC = g*g';
    GCC_A = GCC(A_bin); % vector of GCC scores at connections
    for j = 1:numStats
        gScore(i,j) = corr(q(:,j),GCC_A,'type',whatCorr,'rows','pairwise');
    end
    % Print some info to screen for the user:
    if i==1 || mod(i,numGenes/10)==0
        fprintf(1,'%u/%u\n',i,numGenes);
    end
end

%-------------------------------------------------------------------------------
% Get scores relative to what would be expected from distance
dScores = load('dScores_Spearman.mat','geneEntrez','geneDistanceScores');
[geneEntrez,ia,ib] = intersect(dScores.geneEntrez,[G.GeneStruct.gene_entrez_id]);
gScoresCorrected = gScore(ib) - dScores.geneDistanceScores(ia);

%-------------------------------------------------------------------------------
% Save quickly to .mat file:
save('gScore.mat','gScore','gScoresCorrected','geneEntrez','edgeMeasureNames')

%-------------------------------------------------------------------------------
% Save each result to a separate ErmineJ file:
doCorrected = true;
for i = 1:length(edgeMeasureNames)
    fileName = sprintf('%s',edgeMeasureNames{i});
    if doCorrected
        fileName = [fileName,'_corrected'];
        writeErmineJFile(fileName,gScore(:,i),geneEntrez,edgeMeasureNames{i});
    end
end

% end
