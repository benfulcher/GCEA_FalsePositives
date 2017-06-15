function NullSummaryPlots(GOTable,categoryScores,titleText)

if nargin < 3
    titleText = '';
end
%-------------------------------------------------------------------------------

f = figure('color','w');

% Plot distribution of p-values:
subplot(2,3,1); hold on
histogram(GOTable.pValsZ)
histogram(GOTable.pValsZ_corr)
xlabel('p-values')
legend({'raw','corrected'})
ylabel('frequency')
title(titleText,'interpreter','none')

% Plot distribution of mean nulls:
subplot(2,3,2); hold on
histogram(GOTable.meanNull)
histogram(categoryScores(:,1))
legend('null','real')
plot(ones(2,1)*nanmean(categoryScores(:,1)),[0,max(get(gca,'ylim'))],'r')
xlabel('Category score')
ylabel('Frequency')
title(titleText,'interpreter','none')

% Check dependence on GO category size
subplot(2,3,3)
plot(GOTable.size,GOTable.pValsZ,'.k')
xlabel('GO category size')
ylabel('corrected p-value')
title(titleText,'interpreter','none')

% Relationship between null mean and real scores
subplot(2,3,4); hold on
plot(GOTable.meanNull,categoryScores(:,1),'.k')
plot([min(GOTable.meanNull),max(GOTable.meanNull)],[min(GOTable.meanNull),max(GOTable.meanNull)],'r')
xlabel('mean of null distribution')
ylabel('real scores')
title(titleText,'interpreter','none')

% Std
subplot(2,3,5); hold on
plot(GOTable.stdNull,categoryScores(:,1),'.k')
xlabel('std of null distribution')
ylabel('real scores')
title(titleText,'interpreter','none')

end
