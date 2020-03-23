function RelativeFPSRAutoCorr()
% Does adding spatial autocorrelation to phenotypes boost the correlations to
% GO categories containing genes with high spatial autocorrelation?
%-------------------------------------------------------------------------------

whatSpecies = {'mouse','human'};
whatStructFilt = {'all','cortex'};

GOTable_FPSR = struct();
GOTableCombined = struct();

for s = 1:2
    params = GiveMeDefaultParams(whatSpecies{s},whatStructFilt{s});

    % Get spatial autocorrelation scores as means of individual genes per GO category:
    % (from ComputeSpatialEmbeddingScores)
    load(GiveMeDistanceScoreFileName(params),'GOTable');
    GOTable_ACScores = GOTable;
    clear('GOTable')

    % Get category-level CGE exponential decay parameters
    categoryFileName = sprintf('CategorySpatialScoring_%s-%s.mat',whatSpecies{s},whatStructFilt{s});
    load(categoryFileName,'GOTable')
    GOTable_CategoryLevelScores = GOTable;
    clear('GOTable')

    % Load FPSR (random and spatial autocorrelation)
    % (Ensure that we're not dealing with the reference case):
    params.nulls.customSurrogate = 'none';
    % Load SBP-random FPSR results:
    params.g.whatSurrogate = 'randomMap';
    GOTable_FPSR.(whatSpecies{s}).random = SurrogateEnrichmentProcess(params,false);
    % Load SBP-spatial FPSR results:
    params.g.whatSurrogate = 'spatialLag';
    GOTable_FPSR.(whatSpecies{s}).spatialAC = SurrogateEnrichmentProcess(params,false);

    % Combine/annotate:
    [~,ia,ib] = intersect(GOTable_FPSR.(whatSpecies{s}).random.GOID,GOTable_FPSR.(whatSpecies{s}).spatialAC.GOID);
    GOTableCombined.(whatSpecies{s}) = GOTable_FPSR.(whatSpecies{s}).random(ia,:);
    GOTableCombined.(whatSpecies{s}).sumUnderSigAC = GOTable_FPSR.(whatSpecies{s}).spatialAC.sumUnderSig(ib);

    % Annotate some statistics!:
    % Relative difference from SBP-random -> SBP-spatial:
    CFPR_spatial = GOTableCombined.(whatSpecies{s}).sumUnderSigAC;
    CFPR_rand = GOTableCombined.(whatSpecies{s}).sumUnderSig;
    relDiffFPSR = 100*(CFPR_spatial - CFPR_rand)./CFPR_rand;
    relDiffFPSR_filt = relDiffFPSR;
    relDiffFPSR_filt(CFPR_rand < 10) = nan; % bad statistics
    diffFPSR = CFPR_spatial - CFPR_rand;
    diffFPSR_filt = diffFPSR;
    diffFPSR_filt(CFPR_rand < 10) = nan;
    ratFPSR = CFPR_spatial./CFPR_rand;
    ratFPSR_filt = ratFPSR;
    ratFPSR_filt(CFPR_rand < 10) = nan;
    didIncrease = double(CFPR_spatial > CFPR_rand);
    didIncrease_filt = didIncrease;
    didIncrease_filt(CFPR_rand < 10) = nan;

    GOTableCombined.(whatSpecies{s}).relDiffFPSR = relDiffFPSR;
    GOTableCombined.(whatSpecies{s}).relDiffFPSR_filt = relDiffFPSR_filt;
    GOTableCombined.(whatSpecies{s}).diffFPSR = diffFPSR;
    GOTableCombined.(whatSpecies{s}).diffFPSR_filt = diffFPSR_filt;
    GOTableCombined.(whatSpecies{s}).ratFPSR = ratFPSR;
    GOTableCombined.(whatSpecies{s}).ratFPSR_filt = ratFPSR_filt;
    GOTableCombined.(whatSpecies{s}).didIncrease = didIncrease;
    GOTableCombined.(whatSpecies{s}).didIncrease_filt = didIncrease_filt;

    % Does the level of spatial autocorrelation of genes in a category capture the change?:
    % (annotate the meanScore from GOTable_ACScores)
    [~,ia,ib] = intersect(GOTableCombined.(whatSpecies{s}).GOID,GOTable_ACScores.GOID);
    GOTableCombined.(whatSpecies{s}) = GOTableCombined.(whatSpecies{s})(ia,:);
    GOTableCombined.(whatSpecies{s}).meanACscore = GOTable_ACScores.meanScore(ib);

    % Now let's try to add the category-level scores, GOTable_CategoryLevelScores:
    [~,ia,ib] = intersect(GOTableCombined.(whatSpecies{s}).GOID,GOTable_CategoryLevelScores.GOID);
    GOTableCombined.(whatSpecies{s}) = GOTableCombined.(whatSpecies{s})(ia,:);
    GOTableCombined.(whatSpecies{s}).A_fitted = GOTable_CategoryLevelScores.A_fitted(ib);
    GOTableCombined.(whatSpecies{s}).B_fitted = GOTable_CategoryLevelScores.B_fitted(ib);
    GOTableCombined.(whatSpecies{s}).d0_fitted = GOTable_CategoryLevelScores.d0_fitted(ib);
    GOTableCombined.(whatSpecies{s}).R2fit = GOTable_CategoryLevelScores.R2fit(ib);
    GOTableCombined.(whatSpecies{s}).negRho = GOTable_CategoryLevelScores.negRho(ib);
end

%-------------------------------------------------------------------------------
f = figure('color','w'); hold('on');
ax = gca();
numBins = 10;
theColors = GiveMeColors('mouseHuman');
for s = 1:2
    % isValid = (GOTableCombined.(whatSpecies{s}).R2fit > 0.25);
    % xData = GOTableCombined.(whatSpecies{s}).negRho;
    % xData = GOTableCombined.(whatSpecies{s}).R2fit;
    % xData = GOTableCombined.(whatSpecies{s}).A_fitted; % WORKS WELL FOR MOUSE (NOT HUMAN)
    % xData = GOTableCombined.(whatSpecies{s}).A_fitted - GOTableCombined.(whatSpecies{s}).B_fitted;
    % xData = GOTableCombined.(whatSpecies{s}).d0_fitted;
    xData = GOTableCombined.(whatSpecies{s}).meanACscore;

    % yData = GOTableCombined.(whatSpecies{s}).relDiffFPSR;
    % yData = GOTableCombined.(whatSpecies{s}).relDiffFPSR_filt;
    % ylabel('Relative CFPR(AC) - CFPR(random) (%)')
    % yData = GOTableCombined.(whatSpecies{s}).diffFPSR_filt;
    % yData = GOTableCombined.(whatSpecies{s}).ratFPSR_filt;
    yData = GOTableCombined.(whatSpecies{s}).didIncrease*100;
    ylabel('Increase in CFPR (spatial) (%)')

    BF_PlotQuantiles(xData,yData,numBins,false,false,theColors(s,:),false);
    text(mean(ax.XLim),mean(ax.YLim)*1.1,whatSpecies{s},'Color',theColors(s,:))
    xlabel('Spatial autocorrelation score')

end
f.Position = [1000        1159         239         179];

% Save out as svg file:
fileName = fullfile('OutputPlots','Rel_CFPR_SpatialAC.svg');
saveas(f,fileName,'svg')
fprintf(1,'Saved to %s\n',fileName);

%-------------------------------------------------------------------------------

return % no play for you!
%-------------------------------------------------------------------------------
% PLAYGROUND:
%-------------------------------------------------------------------------------
theSpecies = 'human';
numBins = 10;
isValid = (GOTableCombined.(whatSpecies{s}).R2fit > 0.00);
% BF_PlotQuantiles(GOTableCombined.(theSpecies).d0_fitted(isValid),...
%                     GOTableCombined.(theSpecies).didIncrease(isValid),...
%                     numBins,false,false,theColors(1,:),false);

% [Y,E] = discretize(GOTableCombined.(theSpecies).d0_fitted,numBins);
numDBins = 5;
numThresholds = numDBins+1;
d0Data = GOTableCombined.(theSpecies).d0_fitted(isValid);
xThresholds = arrayfun(@(x)quantile(d0Data,x),linspace(0,1,numThresholds));
xThresholds(end) = xThresholds(end) + eps; % make sure all data included in final bin

f = figure('color','w');
xData = GOTableCombined.(theSpecies).meanACscore(isValid);
yData = GOTableCombined.(theSpecies).didIncrease(isValid)*100;
xLims = zeros(numDBins,2);
yLims = zeros(numDBins,2);
for i = 1:numDBins
    ax = subplot(2,numDBins,i);
    isInBin = (d0Data>=xThresholds(i) & d0Data<xThresholds(i+1));
    BF_PlotQuantiles(xData(isInBin),yData(isInBin),numBins,false,false,theColors(s,:),false);
    title(sprintf('%.1f<d<%.1f',xThresholds(i),xThresholds(i+1)))
    xLims(i,:) = ax.XLim;
    yLims(i,:) = ax.YLim;
end
for i = 1:numDBins
    ax = subplot(2,numDBins,i);
    ax.XLim = [min(xLims(:,1)),max(xLims(:,2))];
    ax.YLim = [min(yLims(:,1)),max(yLims(:,2))];
end
subplot(2,numDBins,numDBins+1:2*numDBins)
histogram(d0Data); %,xThresholds)
title(theSpecies)

% plot(GOTableCombined.(whatSpecies{s}).R2fit,GOTableCombined.(whatSpecies{s}).A_fitted)
%
% histogram(GOTableCombined.(whatSpecies{s}).B_fitted)
%
% scatter(GOTableCombined.mouse.d0_fitted,GOTableCombined.mouse.relDiffFPSR,20,...
%                 GOTableCombined.mouse.meanACscore,'filled')
%
% numBins = 10;
% isInRightDRange = (GOTableCombined.mouse.d0_fitted<3);
% BF_PlotQuantiles(GOTableCombined.mouse.meanACscore(isInRightDRange),...
%                     GOTableCombined.mouse.relDiffFPSR(isInRightDRange),...
%                     numBins,false,false,theColors(1,:),false);
%
% BF_PlotQuantiles(GOTableCombined.mouse.meanACscore(~isInRightDRange),...
%                     GOTableCombined.mouse.relDiffFPSR(~isInRightDRange),...
%                     numBins,false,false,theColors(2,:),false);

end
