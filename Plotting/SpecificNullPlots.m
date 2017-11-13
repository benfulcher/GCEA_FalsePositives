function SpecificNullPlots(GOTable,categoryScores)
%-------------------------------------------------------------------------------
% Look at distribution for some top ones
%-------------------------------------------------------------------------------
numGOCategories = size(categoryScores,1);
numNulls = size(categoryScores,2)-1;
nullInd = 2:numNulls+1;

f = figure('color','w');
for i = 1:min(15,numGOCategories)
    subplot(5,3,i);
    hold on
    histogram(categoryScores(i,nullInd),'edgeColor','k','FaceColor','w');
    plot(categoryScores(i,1)*ones(2,1),[0,max(get(gca,'ylim'))],'-r')
    % plot(whatStat(ix_GO(i))*ones(2,1),[0,max(get(gca,'ylim'))],'-r')
    title(sprintf('%s (%u; p_{corr}=%.2g)\n',GOTable.GOName{i},...
                        GOTable.size(i),GOTable.pValZCorr(i)));
end

end
