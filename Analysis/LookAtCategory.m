
mouseOrHuman = 'mouse';
structFilt = 'cortex';

params = GiveMeDefaultParams(mouseOrHuman,structFilt);
distMat = GiveMeDistanceMatrix(params);

% Get category of interest:
whatGOID = 97369;
[geneData,geneInfo,structInfo,categoryInfo] = GiveMeGOCategory(whatGOID,params);

% Compute 2-d projection of data:
coOrdsXY = mdscale(distMat,2);
% Plot each surrogate spatial map:
f = figure('color','w');
numGenes = height(geneInfo);
for i = 1:numGenes
    subplot(3,5,i)
    mapNorm = zscore(geneData(:,i));
    scatter(coOrdsXY(:,1),coOrdsXY(:,2),25,mapNorm,'filled')
    colormap([flipud(BF_getcmap('blues',9)); 1,1,1; BF_getcmap('reds',9)])
    xlabel('spatialAxis1')
    ylabel('spatialAxis2')
    title(geneInfo.acronym{i})
end

% Where does this category sit in the distribution?:
load(sprintf('CategorySpatialScoring_%s-%s.mat',mouseOrHuman,structFilt),...
                'GOTable');

whatGOID = 18027;

yoDawg = find(GOTable.GOID==whatGOID);

f = figure('color','w');
whatField = 'A_fitted'; %'d0_fitted';
histogram(GOTable.(whatField))
title(sprintf('%s-%f',whatField,GOTable.(whatField)(yoDawg)))

%-------------------------------------------------------------------------------
% Let's look at it:
params.g.whatSurrogate = 'randomMap';
[categoryScores,categoryLabels] = CompareNulls(18027,params,false);
