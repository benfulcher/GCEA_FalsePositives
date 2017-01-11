function [gScore,geneEntrezIDs] = ConnectedGCC(theEdgeData,theGeneData,...
                    geneInfo.entrez_id,corrType,distanceRegressor,absType,...
                    thresholdGoodGene,pValOrStat)
% Idea is to score each gene on it's contribution to a difference between
% connected and unconncted pairs of brain regions
%-------------------------------------------------------------------------------


if nargin < 3
    [GData,geneInfo] = LoadMeG({'robustSigmoid','zscore'},'energy');
    geneEntrezIDs = geneInfo.gene_entrez_id;
end
if nargin < 4
    whatCorr = 'Pearson';
end
if nargin < 5
    distanceRegressor = [];
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
% Get dimensions of gene data
[numRegions,numGenes] = size(geneData);

gScore = zeros(numGenes,1);
parfor i = 1:numGenes
    g = geneData(:,i);
    GCC = g*g';

    % Set to NaN if not enough good values:
    if mean(isnan(GCC(triu(true(size(GCC)),+1)))) > thresholdGoodGene
        gScore(i) = NaN;
        continue
    end

    % Make the two groups:
    GCC_group = cell(2,1);
    if isempty(distanceRegressor)
        % Uncorrected
        GCC_group{1} = GCC(theEdgeData==0); % unconnected
        GCC_group{2} = GCC(theEdgeData==1); % connected
    else
        % Distance regressed out:
        lookyHere = (~isnan(theEdgeData) & ~isnan(GCC));
        [p,S] = polyfit(distanceRegressor(lookyHere),GCC(lookyHere),1);
        GCC_fit = p(2) + p(1)*distanceRegressor;
        GCCresid = GCC - GCC_fit;
    end

    % Do the hypothesis test
    % U-test between connected and unconnected:
    [p,~,stats] = ranksum(GCC_group{1},GCC_group{2});
    % Normalized Mann-Whitney U test (given the sample size may change across features)
    n1 = length(GCC_group{1});
    n2 = length(GCC_group{2});
    normuStat = (stats.ranksum - n1*(n1+1)/2)/n1/n2; % normalized uStat
    switch pValOrStat
    case 'stat'
        gScore(i) = normuStat;
    case 'pVal'
        gScore(i) = pVal;
    end
    if isnan(gScore(i))
        keyboard
    end
end

%-------------------------------------------------------------------------------
% Remove genes with too few good values:
fprintf(1,'Removing %u/%u genes under threshold (too many missing values)\n',...
                    sum(isnan(gScore)),length(gScore));
notGoodEnough = isnan(gScore);
gScore(notGoodEnough) = [];
geneEntrezIDs(notGoodEnough) = [];

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
