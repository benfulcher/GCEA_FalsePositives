function ix_GO = ListCategories(geneInfo,GOTable,numTop,whatStat)

if nargin < 3
    numTop = 30;
end
if nargin < 4
    % Sort GO categories according to this:
    whatStat = 'pValZ'; % meanNull, stdNull, pValZ
end

%-------------------------------------------------------------------------------
numGOCategories = height(GOTable);
numTop = min(numTop,numGOCategories);

[~,ix_GO] = sort(GOTable.(whatStat),'ascend');
fprintf(1,'%u nans removed\n',sum(isnan(GOTable.(whatStat))));
ix_GO(isnan(GOTable.(whatStat)(ix_GO))) = [];

for i = 1:numTop
    geneAcro = geneInfo.acronym(ismember(geneInfo.entrez_id,GOTable.annotations{ix_GO(i)}));
    fprintf(1,'%u (%u genes): %s (p = %.2g; p_corr = %.2g) [%s]\n',...
                i,GOTable.size(ix_GO(i)),...
                GOTable.GOName{ix_GO(i)},...
                GOTable.pValZ(ix_GO(i)),GOTable.pValZCorr(ix_GO(i)),...
                BF_cat(geneAcro));
end

end
