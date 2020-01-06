function RelativeFPSRAutoCorr()
% Does adding spatial autocorrelation to phenotypes boost the correlations to
% GO categories containing genes with high spatial autocorrelation?

numNullSamples_surrogate = 10000;
whatSpecies = {'mouse','human'};

GOTable_FPSR = struct();
GOTableCombined = struct();

for s = 1:2
    params = GiveMeDefaultParams(whatSpecies{s});
    % Compute spatial autocorrelation scores per GO category:
    params.g.normalizationGene = 'zscore';
    params.g.normalizationRegion = 'zscore';
    params.e.numNullSamples = 100; % for speed since we don't actually use the p-values
    results = geneEnrichmentDistance(params);

    % Load FPSR (random and spatial autocorrelation)
    GOTable_FPSR.(whatSpecies{s}).random = SurrogateEnrichmentProcess(whatSpecies{s},numNullSamples_surrogate,'randomUniform','');
    GOTable_FPSR.(whatSpecies{s}).spatialAC = SurrogateEnrichmentProcess(whatSpecies{s},numNullSamples_surrogate,'spatialLag','');

    % Combine/annotate:
    [~,ia,ib] = intersect(GOTable_FPSR.(whatSpecies{s}).random.GOID,GOTable_FPSR.(whatSpecies{s}).spatialAC.GOID);
    GOTableCombined.(whatSpecies{s}) = GOTable_FPSR.(whatSpecies{s}).random(ia,:);
    GOTableCombined.(whatSpecies{s}).sumUnderSigAC = GOTable_FPSR.(whatSpecies{s}).spatialAC.sumUnderSig(ib);
    GOTableCombined.(whatSpecies{s}).relDiffFPSR = 100*(GOTableCombined.(whatSpecies{s}).sumUnderSigAC - ...
                                                        GOTableCombined.(whatSpecies{s}).sumUnderSig)./...
                                                        GOTableCombined.(whatSpecies{s}).sumUnderSig;

    % Does the level of spatial autocorrelation of genes in a category capture the change?:
    % annotate the meanScore
    [~,ia,ib] = intersect(GOTableCombined.(whatSpecies{s}).GOID,results.GOID);
    GOTableCombined.(whatSpecies{s}) = GOTableCombined.(whatSpecies{s})(ia,:);
    GOTableCombined.(whatSpecies{s}).meanACscore = results.meanScore(ib);
end

%-------------------------------------------------------------------------------
f = figure('color','w'); hold('on')
numBins = 10;
theColors = GiveMeColors('mouseHuman');
for s = 1:2
    BF_PlotQuantiles(GOTableCombined.(whatSpecies{s}).meanACscore,GOTableCombined.(whatSpecies{s}).relDiffFPSR,...
                        numBins,false,false,theColors(s,:),false);
    xlabel('Spatial autocorrelation score')
    ylabel('Relative FPSR(AC) - FPSR(random) (%)')
end
f.Position = [1000        1159         239         179];

end
