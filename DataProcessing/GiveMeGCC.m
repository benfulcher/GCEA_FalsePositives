function [gScore,geneEntrezIDs] = GiveMeGCC(edgeData,geneData,geneEntrezIDs,whatCorr,...
                                correctDistance,doAbs,thresholdGoodGene,pValOrStat)
% Returns GCC scores for all genes given some edge metric
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% DEFAULTS:
%-------------------------------------------------------------------------------
if nargin < 3
    [GeneStruct,GData] = LoadMeG(true,{'robustSigmoid','zscore'},'energy');
    geneEntrezIDs = [GeneStruct.gene_entrez_id];
end
if nargin < 4
    whatCorr = 'Pearson';
end
if nargin < 5
    correctDistance = false;
end
if nargin < 6
    doAbs = false;
end
if nargin < 7
    % Need minimum of 50% of data to compute correlation
    thresholdGoodGene = 0.5;
end
if nargin < 8
    pValOrStat = 'stat';
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
    if mean(isnan(GCC_A)) > thresholdGoodGene
        beep
        gScore(i) = NaN;
    else
        switch pValOrStat
        case 'pVal'
            [~,gScore(i)] = corr(edgeVector,GCC_A,'type',whatCorr,'rows','pairwise');
        case 'stat'
            gScore(i) = corr(edgeVector,GCC_A,'type',whatCorr,'rows','pairwise');
        end
    end
    if isnan(gScore(i))
        keyboard
    end
    % Print some info to screen for the user:
    if i==1 || mod(i,numGenes/10)==0
        fprintf(1,'%u/%u\n',i,numGenes);
    end
end

%-------------------------------------------------------------------------------
% Correct for distance:
%-------------------------------------------------------------------------------
% Get scores relative to what would be expected from distance
if correctDistance
    % STILL NEED TO INSPECT WHAT THIS IS DOING...
    if ~strcmp(whatCorr,'Spearman')
        warning('Correcting using Spearman, not %s',whatCorr);
    end
    dScores = load('dScores_Spearman.mat','geneEntrez','geneDistanceScores');
    [geneEntrezMatched,ia,ib] = intersect(dScores.geneEntrez,geneEntrezIDs);
    gScore = bsxfun(@minus,gScore(ib),dScores.geneDistanceScores(ia));
    geneEntrezIDs = geneEntrezMatched;
end

%-------------------------------------------------------------------------------
% Take absolute values:
if doAbs
    gScore = abs(gScore);
end

end
