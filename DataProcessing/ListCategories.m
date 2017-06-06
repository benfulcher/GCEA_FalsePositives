function ListCategories(geneInfo,GOTable,geneEntrezAnnotations,meanNull,pValsZ,pValsZ_corr)


sizeGOCategories = cellfun(@length,geneEntrezAnnotations);
numGOCategories = length(meanNull);

numTop = min(30,numGOCategories);
whatStat = pValsZ; % meanNull, stdNull, pValsZ
[~,ix_GO] = sort(whatStat,'ascend');
fprintf(1,'%u nans removed\n',sum(isnan(whatStat)));
ix_GO(isnan(whatStat(ix_GO))) = [];
for i = 1:numTop
    geneAcro = geneInfo.acronym(ismember(geneInfo.entrez_id,geneEntrezAnnotations{ix_GO(i)}));
    fprintf(1,'%u (%u genes): %s (nullmean = %.2g; p = %.2g; p_corr = %.2g) [%s]\n',i,sizeGOCategories(ix_GO(i)),...
                        GOTable.GOName{ix_GO(i)},meanNull(ix_GO(i)),pValsZ(ix_GO(i)),pValsZ_corr(ix_GO(i)),BF_cat(geneAcro));
end

end
