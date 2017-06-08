function NullSummaryPlots(pValsZ,pValsZ_corr,categoryScores,meanNull,stdNull,sizeGOCategories,titleText)

if nargin < 7
    titleText = '';
end
%-------------------------------------------------------------------------------

f = figure('color','w');

% Plot distribution of p-values:
subplot(2,3,1); hold on
histogram(pValsZ)
histogram(pValsZ_corr)
xlabel('p-values')
legend({'raw','corrected'})
ylabel('frequency')
title(titleText,'interpreter','none')

% Plot distribution of mean nulls:
subplot(2,3,2); hold on
histogram(meanNull)
histogram(categoryScores(:,1))
legend('null','real')
plot(ones(2,1)*nanmean(categoryScores(:,1)),[0,max(get(gca,'ylim'))],'r')
xlabel('Category score')
ylabel('Frequency')
title(titleText,'interpreter','none')

% Check dependence on GO category size
subplot(2,3,3)
plot(sizeGOCategories,pValsZ,'.k')
xlabel('GO category size')
ylabel('corrected p-value')
title(titleText,'interpreter','none')

% Relationship between null mean and real scores
subplot(2,3,4); hold on
plot(meanNull,categoryScores(:,1),'.k')
plot([min(meanNull),max(meanNull)],[min(meanNull),max(meanNull)],'r')
xlabel('mean of null distribution')
ylabel('real scores')
title(titleText,'interpreter','none')

% Std
subplot(2,3,5); hold on
plot(stdNull,categoryScores(:,1),'.k')
xlabel('std of null distribution')
ylabel('real scores')
title(titleText,'interpreter','none')

end
