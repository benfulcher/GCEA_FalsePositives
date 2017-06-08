function SpecificNullPlots()
%-------------------------------------------------------------------------------
% Look at distribution for some top ones
%-------------------------------------------------------------------------------

f = figure('color','w');
for i = 1:min(15,numGOCategories)
    subplot(5,3,i); hold on
    histogram(categoryScores(ix_GO(i),nullInd),'edgeColor','k','FaceColor','w');
    plot(categoryScores(ix_GO(i),1)*ones(2,1),[0,max(get(gca,'ylim'))],'-r')
    % plot(whatStat(ix_GO(i))*ones(2,1),[0,max(get(gca,'ylim'))],'-r')
    title(sprintf('%s (%u; p_{corr}=%.2g)\n',GOTable.GOName{ix_GO(i)},...
                        sizeGOCategories(ix_GO(i)),pValsZ_corr(ix_GO(i))));
    % ,pValsZ(ix(i))
end

end
