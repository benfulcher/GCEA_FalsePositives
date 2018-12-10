function [rawScore,absScore] = IntraCorrelationScore(geneMatrix)

% Compute a pairwise correlation matrix:
C = corr(geneMatrix,'rows','pairwise','type','Pearson');

% Take vector of values in the upper triangle:
corrVect = C(triu(true(size(C)),+1));

% Take the mean of upper triangle as basic statistic:
rawScore = nanmean(corrVect);

% Also take magnitude of any negative correlations (to count negative correlations equally):
absScore = nanmean(abs(corrVect));

end
