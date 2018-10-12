function [Xmean,Xflipped] = flippedMean(X)
% Flips anticorrelated columns to prevent cancellation of signal
%-------------------------------------------------------------------------------

% First find the gene most (abs) correlated to other genes
r = corr(X,'type','Spearman','rows','pairwise');
[~,ix] = max(mean(abs(r)));
clusterCenter = X(:,ix);

% Flip negatively correlated genes
isNeg = (r(ix,:) < 0);
Xflipped = X;
Xflipped(:,isNeg) = -Xflipped(:,isNeg);
fprintf(1,'Flipped %u genes\n',sum(isNeg));

% Then take the mean after flipping
Xmean = nanmean(X,2);

end
