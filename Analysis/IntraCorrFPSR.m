function IntraCorrFPSR()
% Investigate whether FPSR/CFPR relates to intra-category coexpression
%-------------------------------------------------------------------------------

whatSpecies = {'mouse','human'};

% Properties of the intracorrelation statistic:
whatIntraStat = 'raw';
% (these bits about the nulls aren't actually used):
whatShuffle = 'geneShuffle'; % 'geneShuffle', 'independentSpatialShuffle'
numNullSamples_intraCorr = 20000; % (Intra_*_*_20000.mat)

%===============================================================================
% Import (or compute) intra-category correlation data and surrogate random data
%===============================================================================
for s = 1:2
    % Enforce comparison to SBP-random ensemble ('randomMap'):
    params = GiveMeDefaultParams(whatSpecies{s});
    params.g.whatSurrogate = 'randomMap';
    params.nulls.customShuffle = 'none';

    % Load in within-category coexpression stats:
    fileNameIn = sprintf('Intra_%s_%s_%s_%u.mat',whatSpecies{s},whatShuffle,whatIntraStat,numNullSamples_intraCorr);
    resultsIntra = load(fileNameIn);
    fprintf(1,'Imported within-category coexpression information from %s\n',fileNameIn);

    results.(whatSpecies{s}) = struct();
    results.(whatSpecies{s}).intra = resultsIntra.resultsTable;
    results.(whatSpecies{s}).randomReal = SurrogateEnrichmentProcess(params);
end

%-------------------------------------------------------------------------------
% Plot:
%-------------------------------------------------------------------------------
f = figure('color','w'); hold('on');
theXfield = 'intracorr_raw';
theYfield = 'CFPR_random';
numQuantileBins = 10;
theColors = GiveMeColors('mouseHuman');
for s = 1:2
    [~,ia,ib] = intersect(results.(whatSpecies{s}).intra.GOID,results.(whatSpecies{s}).randomReal.GOID);
    GOTableCombined = results.(whatSpecies{s}).intra(ia,:);
    GOTableCombined.CFPR_random = 100*results.(whatSpecies{s}).randomReal.sumUnderSig(ib)/params.nulls.numNullsCFPR;

    BF_PlotQuantiles(GOTableCombined.(theXfield),GOTableCombined.(theYfield),...
            numQuantileBins,false,false,theColors(s,:),false);
end
ylabel('CFPR (SBP-random) (%)')
xlabel(sprintf('Within-category correlation, %s',whatIntraStat))
f.Position = [1000        1159         239         179];

% Save out as svg file:
fileName = fullfile('OutputPlots','IntraCorr_CFPR.svg');
saveas(f,fileName,'svg')
fprintf(1,'Saved to %s\n',fileName);

end
