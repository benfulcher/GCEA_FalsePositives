function GOTable = EstimatePVals(categoryScores,whatTail,GOTable)

if nargin < 3
    whatTail = 'right';
end

%-------------------------------------------------------------------------------
% Summarize nulls, estimate p-values
%-------------------------------------------------------------------------------
numGOCategories = size(categoryScores,1);
nullInd = 2:size(categoryScores,2);
meanNull = nanmean(categoryScores(:,nullInd),2); % mean score of genes in each category
stdNull = nanstd(categoryScores(:,nullInd),[],2); % std of genes in each category
% We should better quantify taking into account the number of nulls:
switch whatTail
case 'right' % categories with higher positive correlations to the edge measure than nulls
    fprintf(1,'Right tail: GO categories with more positive scores than nulls\n');
    pValPerm = arrayfun(@(x)mean(categoryScores(x,nullInd)>=categoryScores(x,1)),...
                                    (1:numGOCategories)');
    pValZ = arrayfun(@(x)1-normcdf(categoryScores(x,1),mean(categoryScores(x,nullInd)),std(categoryScores(x,nullInd))),...
                                    (1:numGOCategories)');
case 'left' % categories with more negative correlations to the edge measure than nulls
    fprintf(1,'Left tail: GO categories with more negative scores than nulls\n');
    pValPerm = arrayfun(@(x)mean(categoryScores(x,nullInd)<=categoryScores(x,1)),...
                                    (1:numGOCategories)');
    pValZ = arrayfun(@(x)normcdf(categoryScores(x,1),mean(categoryScores(x,nullInd)),std(categoryScores(x,nullInd))),...
                                    (1:numGOCategories)');
end

%-------------------------------------------------------------------------------
% Corrected p-values using Benjamini-Hochberg
pValPermCorr = mafdr(pValPerm,'BHFDR',true,'showPlot',false);
pValZCorr = mafdr(pValZ,'BHFDR',true,'showPlot',false);
% q-values of Storey, 2002
% [~,pValZCorr] = mafdr(pValZ);

%-------------------------------------------------------------------------------
% Assign values to categories of GOTable:
%-------------------------------------------------------------------------------
GOTable.pValZ = pValZ;
GOTable.pValZCorr = pValZCorr;
GOTable.pValPerm = pValPerm;
GOTable.pValPermCorr = pValPermCorr;
GOTable.meanScore = categoryScores(:,1);
GOTable.meanNull = meanNull;
GOTable.stdNull = stdNull;

end
