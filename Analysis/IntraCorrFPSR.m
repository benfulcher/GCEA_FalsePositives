function IntraCorrFPSR()
% Investigate whether FPSR relates to intra-category coexpression
%-------------------------------------------------------------------------------

whatSpecies = {'mouse','human'};

% Properties of the intracorrelation statistic:
whatIntraStat = 'raw';
% (these bits about the nulls aren't actually used):
whatShuffle = 'geneShuffle'; % 'geneShuffle', 'independentSpatialShuffle'
numNullSamples_intraCorr = 20000; % (Intra_*_*_20000.mat)

params = GiveMeDefaultParams('mouse');

%===============================================================================
% Import (or compute) intra-category correlation data and surrogate random data
%===============================================================================
for s = 1:2
    results.(whatSpecies{s}) = struct();
    fileNameIn = sprintf('Intra_%s_%s_%s_%u.mat',whatSpecies{s},whatShuffle,whatIntraStat,numNullSamples_intraCorr);
    resultsIntra = load(fileNameIn);
    fprintf(1,'Importing intra-category coexpression data from %s\n',fileNameIn);
    results.(whatSpecies{s}).intra = resultsIntra.resultsTable;
    results.(whatSpecies{s}).randomReal = SurrogateEnrichmentProcess(whatSpecies{s},params.nulls.numNullsFPSR,'randomUniform','');
end

%-------------------------------------------------------------------------------
% Plot:
%-------------------------------------------------------------------------------
f = figure('color','w'); hold('on');
theXfield = 'intracorr_raw';
theYfield = 'FPSR_random';
numQuantileBins = 10;
theColors = GiveMeColors('mouseHuman');
for s = 1:2
    [~,ia,ib] = intersect(results.(whatSpecies{s}).intra.GOID,results.(whatSpecies{s}).randomReal.GOID);
    GOTableCombined = results.(whatSpecies{s}).intra(ia,:);
    GOTableCombined.FPSR_random = 100*results.(whatSpecies{s}).randomReal.sumUnderSig(ib)/params.nulls.numNullsFPSR;

    BF_PlotQuantiles(GOTableCombined.(theXfield),GOTableCombined.(theYfield),...
            numQuantileBins,false,false,theColors(s,:),false);
end
ylabel('FPSR-random (%)')
xlabel(sprintf('Intra-category correlation, %s',whatIntraStat))
f.Position = [1000        1159         239         179];

end
