function gScore = GiveMeGCC(edgeData,geneData,whatCorr,correctDistance,doAbs,thresholdGoodGene)
% Returns GCC scores for all genes given some edge metric
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% DEFAULTS:
%-------------------------------------------------------------------------------
if nargin < 2
    [~,GData] = LoadMeG(true,{'robustSigmoid','zscore'},'energy');
end
if nargin < 3
    whatCorr = 'Pearson';
end
if nargin < 4
    correctDistance = false;
end
if nargin < 5
    doAbs = false;
end
if nargin < 6
    % Need minimum of 50% of data to compute correlation
    thresholdGoodGene = 0.5;
end
%-------------------------------------------------------------------------------

[numRegions,numGenes] = size(geneData);

% Convert to vector across edges:
isEdge = (edgeData~=0);
fprintf(1,'%u edges in the data\n',sum(isEdge(:)));
edgeVector = edgeData(isEdge);

%-------------------------------------------------------------------------------
% Ok, so now we can find correlations to GCC scores across genes
gScore = zeros(numGenes,1);
fprintf(1,'Looping over %u genes\n',numGenes);
for i = 1:numGenes
    g = geneData(:,i);
    GCC = g*g';
    GCC_A = GCC(isEdge);
    if mean(isnan(GCC_A)) < thresholdGoodGene
        gScore(i) = NaN;
    else
        gScore(i) = corr(edgeVector,GCC_A,'type',whatCorr,'rows','pairwise');
    end
    % Print some info to screen for the user:
    if i==1 || mod(i,numGenes/10)==0
        fprintf(1,'%u/%u\n',i,numGenes);
    end
end

%-------------------------------------------------------------------------------
% Take absolute values:
if doAbs
    gScore = abs(gScore);
end

%-------------------------------------------------------------------------------
% Correct for distance:
%-------------------------------------------------------------------------------
% Get scores relative to what would be expected from distance
dScores = load('dScores_Spearman.mat','geneEntrez','geneDistanceScores');
[geneEntrezMatched,ia,ib] = intersect(dScores.geneEntrez,geneEntrezIDs);
gScoresCorrected = bsxfun(@minus,gScore(ib,:),dScores.geneDistanceScores(ia));

end
