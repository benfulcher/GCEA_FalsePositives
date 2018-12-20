% CorrCoexpNullWidth
%-------------------------------------------------------------------------------

numGOIDs = 50;
whatSpecies = 'mouse';
whatSurrogate = 'spatialShuffle';
% whatSurrogate = 'spatialLag';
numNullSamplesIntraEstimate = 100;
numNullSamplesNullWidth = 500;
whatCorr = 'Spearman';
theCategorySize = 20; % restrict GO categories to being this size.

%-------------------------------------------------------------------------------
% Get default parameters:
mouseParams = GiveMeDefaultParams(whatSpecies);


%-------------------------------------------------------------------------------
% Pick GO IDs across a range of mean intra-category coexpression
% (need to use significance testing by category size)


GO_intraCoexp = IntraCorrelationByCategory(params,'geneShuffle',numNullSamplesIntraEstimate);
theRange = round(linspace(1,height(GO_intraCoexp),numGOIDs));
GOIDs = GO_intraCoexp.GOID(theRange);
GONames = GO_intraCoexp.GOName(theRange);
pValZCorrs = GO_intraCoexp.pValZCorr(theRange);

%-------------------------------------------------------------------------------
% Estimate null width under the specified null model
params = GiveMeDefaultParams(whatSpecies);
params.g.whatSurrogate = whatSurrogate;

nullWidth = zeros(numGOIDs);
for i = 1:numGOIDs
    categoryScores = GiveMeCategoryNullDist(GOIDs(i),params,numNullSamplesNullWidth,whatCorr);
    nullWidth(i) = var(categoryScores);
end

%-------------------------------------------------------------------------------
% Plot:
f = figure('color','w');
plot(pValZCorrs,nullWidth,'ok');
