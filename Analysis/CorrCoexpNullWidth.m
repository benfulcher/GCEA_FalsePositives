% CorrCoexpNullWidth
%-------------------------------------------------------------------------------

numGOIDs = 50;
whatSpecies = 'mouse';
whatSurrogate = 'spatialShuffle';
% whatSurrogate = 'spatialLag';
numNullSamplesIntraEstimate = 100;
numNullSamplesNullWidth = 500;
whatCorr = 'Spearman';
minAnnotations = 20; % restrict GO categories to being this size.

%-------------------------------------------------------------------------------
% Get default parameters:
params = GiveMeDefaultParams(whatSpecies);
% only consider categories with at least theCategorySize annotations
params.e.sizeFilter = [minAnnotations,200];
params.e.sizeFix = minAnnotations;

%-------------------------------------------------------------------------------
% Pick GO IDs across a range of mean intra-category coexpression
% (need to use significance testing by category size)
GO_intraCoexp = IntraCorrelationByCategory(params,'geneShuffle',numNullSamplesIntraEstimate);
theRange = round(linspace(1,height(GO_intraCoexp),numGOIDs));
GOIDs = GO_intraCoexp.GOID(theRange);
GONames = GO_intraCoexp.GOName(theRange);
pValZCorrs = GO_intraCoexp.pValZCorr(theRange);
meanIntraCorr = GO_intraCoexp.mouse(theRange);

%-------------------------------------------------------------------------------
% Estimate null width under the specified null model
params.g.whatSurrogate = whatSurrogate;

nullWidth = zeros(numGOIDs,1);
for i = 1:numGOIDs
    fprintf(1,'\n\nNull %u/%u\n\n',i,numGOIDs);
    categoryScores = GiveMeCategoryNullDist(GOIDs(i),params,numNullSamplesNullWidth,whatCorr);
    nullWidth(i) = var(categoryScores);
end

%-------------------------------------------------------------------------------
% Plot:
f = figure('color','w');
plot(meanIntraCorr,nullWidth,'ok');
xlabel('mean intra-category correlation')
ylabel(sprintf('%s:nullWidth',whatSurrogate))
[r,p] = corr(meanIntraCorr,nullWidth);
title(sprintf('r = %.3g, p = %.3g',r,p))
