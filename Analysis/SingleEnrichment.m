function GOTable = SingleEnrichment(geneEntrezIDs,geneScores,processFilter,sizeFilter,numIters)
%-------------------------------------------------------------------------------
% Do an ermineJ style analysis for a given set of entrezIDs and scores
%-------------------------------------------------------------------------------

if nargin < 3
    processFilter = 'biological_process';
end
if nargin < 4
    sizeFilter = [5,200];
end
if nargin < 5
    numIters = 10000;
end
%-------------------------------------------------------------------------------
numGenes = length(geneScores);

% Retrieve GO annotations:
[GOTable,geneEntrezAnnotations] = GetFilteredGOData(processFilter,sizeFilter,geneEntrezIDs);
sizeGOCategories = cellfun(@length,geneEntrezAnnotations);
% Add to table:
GOTable.size = sizeGOCategories;
numGOCategories = height(GOTable);

%-------------------------------------------------------------------------------
% Compute the mean score for within each category:
categoryScores = zeros(numGOCategories,1);
for j = 1:numGOCategories
    matchMe = ismember(geneEntrezIDs,geneEntrezAnnotations{j});
    if sum(matchMe) <= 1
        continue
    end
    categoryScores(j) = nanmean(geneScores(matchMe));
end

%-------------------------------------------------------------------------------
% Generate a null distribution
uniqueSizes = unique(sizeGOCategories);
numSizes = length(uniqueSizes);
nullDistribution = zeros(numSizes,numIters);
for j = 1:numIters
    rp = randperm(numGenes); % takes a millisecond
    for i = 1:numSizes
        try
            nullDistribution(i,j) = nanmean(geneScores(rp(1:uniqueSizes(i))));
        catch
            keyboard
        end
    end
end

%-------------------------------------------------------------------------------
% Compute p-values
pVals = zeros(numGOCategories);
for i = 1:numGOCategories
    pVals(i) = mean(categoryScores(i) > nullDistribution(uniqueSizes==sizeGOCategories(i),:));
end

%-------------------------------------------------------------------------------
% FDR correct:
pVals_corr = mafdr(pVals,'BHFDR','true');

%-------------------------------------------------------------------------------
% Update the GO table:
GOTable.pVal = pVals;
GOTable.pVal_corr = pVals_corr;

end
