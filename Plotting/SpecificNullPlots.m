function SpecificNullPlots(categoryScores,GOTable,ix_GO)
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
    histogram(categoryScores(ix_GO(i),nullInd),'edgeColor','k','FaceColor','w');
    plot(categoryScores(ix_GO(i),1)*ones(2,1),[0,max(get(gca,'ylim'))],'-r')
    % plot(whatStat(ix_GO(i))*ones(2,1),[0,max(get(gca,'ylim'))],'-r')
    title(sprintf('%s (%u; p_{corr}=%.2g)\n',GOTable.GOName{ix_GO(i)},...
                        GOTable.size(ix_GO(i)),GOTable.pValZ_corr(ix_GO(i))));
end

end
