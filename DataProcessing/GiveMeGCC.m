function [gScore,geneEntrezIDs] = GiveMeGCC(edgeData,geneData,geneEntrezIDs,whatCorr,...
                                correctDistance,absType,thresholdGoodGene,pValOrStat)
% Returns GCC scores for all genes given some edge metric
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% DEFAULTS:
%-------------------------------------------------------------------------------
if nargin < 3
    [GData,geneInfo] = LoadMeG({'robustSigmoid','zscore'},'energy');
    geneEntrezIDs = geneInfo.gene_entrez_id;
end
if nargin < 4
    whatCorr = 'Pearson';
end
if nargin < 5
    correctDistance = false;
end
if nargin < 6
    absType = 'pos'; % measure how GCC increases with the edge statistic
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
if any(size(edgeData)==1)
    error('edgeData must be a square matrix');
end
isEdge = (edgeData~=0);
fprintf(1,'%u edges in the data\n',sum(isEdge(:)));
edgeVector = edgeData(isEdge);

%-------------------------------------------------------------------------------
% Ok, so now we can find correlations to GCC scores across genes
gScore = zeros(numGenes,1);
fprintf(1,'Looping over %u genes, computing correlations across %u edges\n',...
                                                numGenes,length(edgeVector));
parfor i = 1:numGenes
    g = geneData(:,i);
    GCC = g*g';
    GCC_A = GCC(isEdge);
    if mean(isnan(GCC_A)) > thresholdGoodGene
        gScore(i) = NaN;
    else
        % Compute the correlation statistic:
        [rho,pVal] = corr(edgeVector,GCC_A,'type',whatCorr,'rows','pairwise');
        switch pValOrStat
        case 'pVal'
            gScore(i) = pVal;
        case 'stat'
            gScore(i) = rho;
        end
        if isnan(gScore(i))
            keyboard
        end
    end
    % Print some info to screen for the user [not useful with parfor]:
    % if i==1 || mod(i,round(numGenes/5))==0
    %     fprintf(1,'%u/%u\n',i,numGenes);
    % end
end

% Remove genes with too few good values:
fprintf(1,'Removing %u/%u genes under threshold (too many missing values)\n',...
                    sum(isnan(gScore)),length(gScore));
notGoodEnough = isnan(gScore);
gScore(notGoodEnough) = [];
geneEntrezIDs(notGoodEnough) = [];

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
% Transform the scores:
switch absType
case 'pos'
    fprintf(1,'Scores capture how GCC *increases* with the edge statistic\n');
case 'neg'
    gScore = -gScore;
    fprintf(1,'Scores converted to negatives (measuring how much GCC *decreases* with edge statistic)\n');
case 'abs'
    gScore = abs(gScore);
    fprintf(1,'Scores converted to absolute values -- measure how GCC correlates with the edge statistic (+ly or -ly equivalent)\n');
otherwise
    error('Unknown abs type: %s',absType);
end

end
