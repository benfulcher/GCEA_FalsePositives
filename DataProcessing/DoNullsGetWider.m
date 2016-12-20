% DoNullsGetWider
%-------------------------------------------------------------------------------
% Idea is to check intuition that GO categories containing genes that
% are spatially correlated across the brain should have nulls that are
% wider, as genes move together when there is some random fluctuation
% so they should on average move further from zero
%-------------------------------------------------------------------------------
whatNull = 'uniformTopology'; %'permutedGeneDep';


corrs = load('whatGenesCorrelate-Spearman-biological_process-Gnone_Rnone.mat');
meanUpper = @(corrMat) mean(corrMat(triu(true(size(corrMat)),+1)));
meanCorrs = cellfun(meanUpper,corrs.PC);

switch whatNull
case 'permutedGeneDep'
    nulls = load('ktotktot-permutedGeneDep-biological_process-Gnone_Rnone-250nulls.mat');
case 'uniformTopology'
    nulls = load('ktotktot-biological_process-250nulls.mat');
end
nullInd = 2:nulls.numNulls+1;
stdNull = nanstd(nulls.categoryScores(:,nullInd),[],2); % std of genes in each category
meanNull = nanmean(nulls.categoryScores(:,nullInd),2);

%-------------------------------------------------------------------------------
% Match tables:
%-------------------------------------------------------------------------------
ids1 = corrs.GOTable.GOID;
ids2 = nulls.GOTable.GOID;
[~,ia,ib] = intersect(ids1,ids2);
meanCorrs = meanCorrs(ia);
meanNull = meanNull(ib);
stdNull = stdNull(ib);

%-------------------------------------------------------------------------------
% Plot relationship
%-------------------------------------------------------------------------------
f = figure('color','w');
subplot(121)
plot(meanCorrs,stdNull,'.k');
r = corr(meanCorrs,stdNull);
xlabel('mean pairwise correlation within gene group')
ylabel('variance of dependent gene permuted nulls')
title(sprintf('%s: r = %.2f',whatNull,r))
subplot(122)
plot(meanCorrs,meanNull,'.k');
r = corr(meanCorrs,meanNull);
xlabel('mean pairwise correlation within gene group')
ylabel('mean of dependent gene permuted nulls')
title(sprintf('%s: r = %.2f',whatNull,r))
