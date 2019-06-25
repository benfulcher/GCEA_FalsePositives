function [rawScore,absScore,VE1] = IntraCorrelationScore(geneMatrix)

% Compute a pairwise correlation matrix:
C = corr(geneMatrix,'rows','pairwise','type','Pearson');

% Take vector of values in the upper triangle:
corrVect = C(triu(true(size(C)),+1));

% Take the mean of upper triangle as basic statistic:
rawScore = nanmean(corrVect);

% Also take magnitude of any negative correlations (to count negative correlations equally):
absScore = nanmean(abs(corrVect));

%-------------------------------------------------------------------------------
% Alternative score is the variance explained by the first principal component:
normGeneMatrix = BF_NormalizeMatrix(geneMatrix,'zscore');
if any(isnan(normGeneMatrix))
    [coeff,score,latent] = pca(normGeneMatrix,'Algorithm','als');
else
    [coeff,score,latent] = pca(normGeneMatrix);
end
perc = latent/sum(latent)*100;
VE1 = perc(1);
if isnan(VE1)
    keyboard
end

end
