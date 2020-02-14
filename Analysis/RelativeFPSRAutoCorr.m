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
    GOTableCombined.(whatSpecies{s}).relDiffFPSR = 100*(GOTableCombined.(whatSpecies{s}).sumUnderSigAC - ...
                                                        GOTableCombined.(whatSpecies{s}).sumUnderSig)./...
                                                        GOTableCombined.(whatSpecies{s}).sumUnderSig;

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
f = figure('color','w'); hold('on')
numBins = 10;
theColors = GiveMeColors('mouseHuman');
for s = 1:2
    % isValid = (GOTableCombined.(whatSpecies{s}).R2fit > 0.2);
    % xData = GOTableCombined.(whatSpecies{s}).A_fitted-GOTableCombined.(whatSpecies{s}).B_fitted;
    xData = GOTableCombined.(whatSpecies{s}).meanACscore;
    BF_PlotQuantiles(xData,...
                        GOTableCombined.(whatSpecies{s}).relDiffFPSR,...
                        numBins,false,false,theColors(s,:),false);
    xlabel('Spatial autocorrelation score')
    ylabel('Relative FPSR(AC) - FPSR(random) (%)')
end
f.Position = [1000        1159         239         179];

%-------------------------------------------------------------------------------
% PLAYGROUND:
%-------------------------------------------------------------------------------
% BF_PlotQuantiles(GOTableCombined.mouse.d0_fitted,...
%                     GOTableCombined.mouse.relDiffFPSR,...
%                     numBins,false,false,theColors(1,:),false);
%
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
