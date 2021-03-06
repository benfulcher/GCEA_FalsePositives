function WhatGeneGroupsCorrelate()

% Params:
whatCorr = 'Spearman';
thresholdGoodGene = 0.5; % threshold of valid coexpression values at which a gene is kept

% Get gene data:
energyOrDensity = 'energy'; % what gene expression data to use
normalizationGene = 'none'; % 'none', 'mixedSigmoid'
normalizationRegion = 'none'; % 'none', 'zscore'
[geneData,geneInfo,structInfo] = LoadMeG({normalizationGene,normalizationRegion},energyOrDensity);

processFilter = 'biological_process';
sizeFilter = [5,200];
GOTable = GetFilteredGOData(processFilter,sizeFilter,geneInfo.entrez_id);
numGOCategories = height(GOTable);

%-------------------------------------------------------------------------------
% Go through and compute pairwise correlations within each category:
PC = cell(numGOCategories,1);
parfor i = 1:numGOCategories
    matchMe = ismember(geneInfo.entrez_id,GOTable.annotations{i});
    if sum(matchMe) <= 1;
        continue
    end
    % Compute pairwise correlation of expression for genes in this category:
    geneSubset = geneData(:,matchMe);
    goodEnough = mean(isnan(geneSubset)) < thresholdGoodGene;
    if ~all(goodEnough)
        geneSubset = geneSubset(:,goodEnough);
    end
    PC{i} = corr(geneSubset,'rows','pairwise','type',whatCorr);

    if i==1 || mod(i,round(numGOCategories/10))==0
        fprintf(1,'%u/%u\n',i,numGOCategories);
    end
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Save
%-------------------------------------------------------------------------------
fileName = sprintf('whatGenesCorrelate-%s-%s-G%s_R%s.mat',whatCorr,...
                processFilter,normalizationGene,normalizationRegion);
save(fullfile('DataOutputs',fileName));
fprintf(1,'Saved %s\n',fileName);

%-------------------------------------------------------------------------------
% Find mean correlation
meanUpper = @(corrMat) mean(corrMat(triu(true(size(corrMat)),+1)));
meanCorrs = cellfun(meanUpper,PC);

% What are the top categories:
[~,ix] = sort(meanCorrs,'descend');
ix(isnan(meanCorrs(ix))) = [];
% toNumber = @(GOCell) cellfun(@(x)str2num(x(4:end)),GOCell,'UniformOutput',true);
% tableIDs = arrayfun(@(x)str2num(x(4:end)),annotationTable.GO);

numTopCategories = 100;
for i = 1:numTopCategories
    % Match to table:
    fprintf(1,'%u (%u genes): %s (%.2f)\n',i,GOTable.size(ix(i)),GOTable.GOName{ix(i)},meanCorrs(ix(i)));
end

%-------------------------------------------------------------------------------
% Explore:
f = figure('color','w');
plot(sizeGOCategories,meanCorrs,'.k');
xlabel('Category size')
ylabel(sprintf('mean %s correlation of genes in category',whatCorr))

end
