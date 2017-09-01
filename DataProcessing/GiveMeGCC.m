function [gScore,geneEntrezIDs] = GiveMeGCC(edgeData,geneData,geneEntrezIDs,whatCorr,...
                                distanceRegressor,absType,thresholdGoodGene,pValOrStat)
% Returns GCC scores for all genes given some edge metric
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% DEFAULTS:
%-------------------------------------------------------------------------------
if nargin < 3
    [geneData,geneInfo] = LoadMeG({'robustSigmoid','zscore'},'energy');
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

%-------------------------------------------------------------------------------
% Ok, so now we can find correlations to GCC scores across genes
gScore = zeros(numGenes,1);
if strcmp(whatCorr,'ttest') % Group based on edge data = {0,1}
    % Some info for the user:
    if isempty(distanceRegressor)
        fprintf(1,'***COMPUTING CORRELATIONS (WITHOUT A DISTANCE REGRESSOR)^^^\n');
    else
        fprintf(1,'***COMPUTING PARTIAL CORRELATIONS (WITH A LINEAR LEAST SQUARES DISTANCE REGRESSOR)***\n');
    end
    fprintf(1,'~~~%u (group 1) versus %u (group 2)\n',sum(edgeData(:)==1),sum(edgeData(:)==0));

    % Loop over genes to compute scores:
    for i = 1:numGenes
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
            GCC_group{1} = GCC(edgeData==1 & ~isnan(GCC)); % connected
            GCC_group{2} = GCC(edgeData==0 & ~isnan(GCC)); % unconnected
        else
            % Distance regressed out of GCC scores:
            lookyHere = (~isnan(edgeData) & ~isnan(GCC));
            p = polyfit(distanceRegressor(lookyHere),GCC(lookyHere),1);
            GCC_fit = p(2) + p(1)*distanceRegressor; % linear fit GCC to distance
            GCCresid = GCC - GCC_fit; % residuals from linear fit
            GCC_group{1} = GCCresid(edgeData==1 & ~isnan(GCCresid)); % connected
            GCC_group{2} = GCCresid(edgeData==0 & ~isnan(GCCresid)); % unconnected
        end

        % Compute the gene score
        gScore(i) = (mean(GCC_group{1})-mean(GCC_group{2}))/sqrt(std(GCC_group{1})^2 + std(GCC_group{2})^2); %/...
                % (sqrt((std(GCC_group{1})/length(GCC_group{1}))+sqrt((std(GCC_group{2})/length(GCC_group{2})))));
        % [h,pVal,~,stats] = ttest2(GCC_group{1},GCC_group{2},'Vartype','unequal');
        % % % U-test between connected and unconnected:
        % [p,~,stats] = ranksum(GCC_group{1},GCC_group{2});
        % % Normalized Mann-Whitney U test (given the sample size may change across features)
        % n1 = length(GCC_group{1});
        % n2 = length(GCC_group{2});
        % normuStat = (stats.ranksum - n1*(n1+1)/2)/n1/n2; % normalized uStat
        % switch pValOrStat
        % case 'stat'
        %     gScore(i) = stats.tstat;
        %     % gScore(i) = normuStat;
        % case 'pVal'
        %     gScore(i) = pVal;
        % end
        % if isnan(gScore(i))
        %     keyboard
        % end
    end
else
    % Correlations with the outer product and the edge vector computed across edges
    if any(size(edgeData)==1)
        error('edgeData must be a square matrix');
    end

    isEdge = (edgeData~=0);
    % Convert to vector across edges:
    fprintf(1,'Computing scores only across %u edges in the data\n',sum(isEdge(:)));
    edgeVector = edgeData(isEdge);
    if ~isempty(distanceRegressor)
        fprintf(1,'***COMPUTING PARTIAL CORRELATIONS USING DISTANCE AS A REGRESSOR***\n');
        distanceRegressor = distanceRegressor(isEdge);
    else
        fprintf(1,'***COMPUTING CORRELATIONS (WITHOUT ANY REGRESSORS)^^^\n');
    end


    fprintf(1,'Looping over %u genes, computing %s correlations across %u edges...\n',...
                                                numGenes,whatCorr,length(edgeVector));
    parfor i = 1:numGenes
        g = geneData(:,i);
        GCC = g*g';
        GCCVector = GCC(isEdge);
        if mean(isnan(GCCVector)) > thresholdGoodGene
            gScore(i) = NaN;
        else
            if ~isempty(distanceRegressor)
                % Partial correlations using distance as a regressor
                [rho,pVal] = partialcorr([edgeVector,GCCVector],distanceRegressor,...
                                    'rows','pairwise','type',whatCorr);
                switch pValOrStat
                case 'pVal'
                    gScore(i) = pVal(1,2);
                case 'stat'
                    gScore(i) = rho(1,2);
                end
            else
                % Pure correlations:
                [rho,pVal] = corr(edgeVector,GCCVector,'type',whatCorr,'rows','pairwise');
                switch pValOrStat
                case 'pVal'
                    gScore(i) = pVal;
                case 'stat'
                    gScore(i) = rho;
                end
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
end

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


%-------------------------------------------------------------------------------
% Correct for distance:
%-------------------------------------------------------------------------------
% Get scores relative to what would be expected from distance
% A pretty poor way to do it because the signature of distance may be from a
% different pattern of GCC scores than the edge metric of interest
% if correctDistance
%     % STILL NEED TO INSPECT WHAT THIS IS DOING...
%     if ~strcmp(whatCorr,'Spearman')
%         warning('Correcting using Spearman, not %s',whatCorr);
%     end
%     dScores = load('dScores_Spearman.mat','geneEntrez','geneDistanceScores');
%     [geneEntrezMatched,ia,ib] = intersect(dScores.geneEntrez,geneEntrezIDs);
%     gScore = bsxfun(@minus,gScore(ib),dScores.geneDistanceScores(ia));
%     geneEntrezIDs = geneEntrezMatched;
% end


% Plot some of the top ones:
% geneData(:,notGoodEnough) = [];
% [~,ix] = sort(gScore,'descend');
% f = figure('color','w');
% for i = 1:10
%     subplot(2,5,i);
%     g = geneData(:,ix(i));
%     GCC = g*g';
%     GCCVector = GCC(isEdge);
%     plot(edgeVector,GCCVector,'.k');
%     title(sprintf('%.2f -- %u',gScore(ix(i)),geneEntrezIDs(ix(i))));
%     xlabel('edge properties')
%     ylabel('gcc')
% end
%
% % Look at distributions of genes that did well:
% f = figure('color','w');
% for i = 1:20
%     subplot(5,8,i)
%     histogram(geneData(:,ix(i)))
%     title(i)
% end
% for i = 1:20
%     subplot(5,8,20+i)
%     histogram(geneData(:,ix(end-i)))
%     title(length(ix)-i)
% end
% keyboard


end
