function [gScore,geneEntrezIDs] = GiveMeGCC(edgeData,geneData,geneEntrezIDs,params,dRegressor)
% Returns GCC scores for all genes given some edge metric
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs & set defaults:
if nargin < 3
    [geneData,geneInfo] = LoadMeG(params.g)
    geneEntrezIDs = geneInfo.gene_entrez_id;
end
if nargin < 4
    params = GiveMeDefaultParams('mouse');
end
if nargin < 5
    dRegressor = []; % don't use a distance regressor
end

%-------------------------------------------------------------------------------
% Get dimensions of gene data
[numRegions,numGenes] = size(geneData);

%-------------------------------------------------------------------------------
% Check expression data ranges?
if any(geneData(:) < 0)
    warning('GCC scores don''t make too much sense when data have negatives... :-O')
end

%-------------------------------------------------------------------------------
% Ok, so now we can find correlations to GCC scores across genes
gScore = zeros(numGenes,1);
if strcmp(params.gcc.whatCorr,'ttest')
    % Group based on edge data = {0,1}
    if isempty(dRegressor)
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
        if mean(isnan(GCC(triu(true(size(GCC)),+1)))) > params.gcc.thresholdGoodGene
            gScore(i) = NaN;
            continue
        end

        % Make the two groups:
        GCC_group = cell(2,1);
        if isempty(dRegressor)
            % Uncorrected
            GCC_group{1} = GCC(edgeData==1 & ~isnan(GCC)); % connected
            GCC_group{2} = GCC(edgeData==0 & ~isnan(GCC)); % unconnected
        else
            % Distance regressed out of GCC scores:
            lookyHere = (~isnan(edgeData) & ~isnan(GCC));
            p = polyfit(dRegressor(lookyHere),GCC(lookyHere),1);
            GCC_fit = p(2) + p(1)*dRegressor; % linear fit GCC to distance
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
    if ~isempty(dRegressor)
        fprintf(1,'***COMPUTING PARTIAL CORRELATIONS USING DISTANCE AS A REGRESSOR***\n');
        dRegressor = dRegressor(isEdge);
    else
        fprintf(1,'***COMPUTING CORRELATIONS (WITHOUT ANY REGRESSORS)^^^\n');
    end


    fprintf(1,'Looping over %u genes, computing %s correlations across %u edges...\n',...
                                                numGenes,params.gcc.whatCorr,length(edgeVector));
    parfor i = 1:numGenes
        g = geneData(:,i);
        GCC = g*g'; % self product at each edge
        GCCVector = GCC(isEdge);
        if mean(isnan(GCCVector)) > params.gcc.thresholdGoodGene
            gScore(i) = NaN;
        else
            if isempty(dRegressor)
                % Pure correlations:
                [rho,pVal] = corr(edgeVector,GCCVector,'type',params.gcc.whatCorr,...
                                                'rows','pairwise');
                switch params.gcc.pValOrStat
                case 'pVal'
                    gScore(i) = pVal;
                case 'stat'
                    gScore(i) = rho;
                end
            else
                % Partial correlations using distance as a regressor
                [rho,pVal] = partialcorr([edgeVector,GCCVector],dRegressor,...
                                    'rows','pairwise','type',params.gcc.whatCorr);
                switch params.gcc.pValOrStat
                case 'pVal'
                    gScore(i) = pVal(1,2);
                case 'stat'
                    gScore(i) = rho(1,2);
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
if any(isnan(gScore))
    fprintf(1,'Removing %u/%u genes under threshold (too many missing values)\n',...
                        sum(isnan(gScore)),length(gScore));
    notGoodEnough = isnan(gScore);
    gScore(notGoodEnough) = [];
    geneEntrezIDs(notGoodEnough) = [];
else
    fprintf(1,'All genes had good enough values for this analysis\n');
end

%-------------------------------------------------------------------------------
% Transform the scores:
fprintf(1,'Across %u genes, mean score is %.3f\n',length(gScore),mean(gScore));
switch params.gcc.absType
case 'pos'
    fprintf(1,'Scores capture how GCC *increases* with the edge statistic\n');
case 'neg'
    gScore = -gScore;
    fprintf(1,'Scores converted to negatives (measuring how much GCC *decreases* with edge statistic)\n');
case 'abs'
    gScore = abs(gScore);
    fprintf(1,'Scores converted to absolute values -- measure how GCC correlates with the edge statistic (+ly or -ly equivalent)\n');
otherwise
    error('Unknown abs type: %s',params.gcc.absType);
end
fprintf(1,'After transformation, mean score is %.3f\n',mean(gScore));


end
