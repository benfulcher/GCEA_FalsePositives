function ix_GO = ListCategories(geneInfo,GOTable,geneEntrezAnnotations,numTop)

if nargin < 7
    numTop = 30;
end

% Sort GO categories according to this:
whatStat = 'pValsZ'; % meanNull, stdNull, pValsZ

%-------------------------------------------------------------------------------
numGOCategories = height(GOTable);
numTop = min(numTop,numGOCategories);

[~,ix_GO] = sort(GOTable.(whatStat),'ascend');
fprintf(1,'%u nans removed\n',sum(isnan(GOTable.(whatStat))));
ix_GO(isnan(GOTable.(whatStat)(ix_GO))) = [];

for i = 1:numTop
    geneAcro = geneInfo.acronym(ismember(geneInfo.entrez_id,geneEntrezAnnotations{ix_GO(i)}));
    fprintf(1,'%u (%u genes): %s (nullmean = %.2g; p = %.2g; p_corr = %.2g) [%s]\n',...
                i,GOTable.size(ix_GO(i)),...
                GOTable.GOName{ix_GO(i)},GOTable.meanNull(ix_GO(i)),...
                GOTable.pValsZ(ix_GO(i)),GOTable.pValsZ_corr(ix_GO(i)),...
                BF_cat(geneAcro));
end

end
