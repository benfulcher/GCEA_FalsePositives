function IntraCorrFPSR()

whatSpecies = {'mouse','human'};
whatShuffle = 'geneShuffle'; % 'geneShuffle', 'independentSpatialShuffle'
whatIntraStat = 'raw';
numNullSamples_intraCorr = 20000; % (Intra_*_*_20000.mat)
numNullSamples_surrogate = 10000; % (SurrogateGOTables_10000_*.mat)

%===============================================================================
% Import (or compute) intra-category correlation data and surrogate random data
%===============================================================================
for s = 1:2
    results.(whatSpecies{s}) = struct();
    fileNameIn = sprintf('Intra_%s_%s_%s_%u.mat',whatSpecies{s},whatShuffle,whatIntraStat,numNullSamples_intraCorr);
    resultsIntra = load(fileNameIn);
    fprintf(1,'Importing intra-category coexpression data from %s\n',fileNameIn);
    results.(whatSpecies{s}).intra = resultsIntra.resultsTable;
    results.(whatSpecies{s}).randomReal = SurrogateEnrichmentProcess(whatSpecies{s},numNullSamples_surrogate,'randomUniform','');
end

%-------------------------------------------------------------------------------
% Plot:
%-------------------------------------------------------------------------------
f = figure('color','w'); hold('on');
theXfield = 'intracorr_raw';
theYfield = 'FPSR_random';
numQuantileBins = 10;
spectralColors = BF_getcmap('spectral',3,1);
theColors = {'k',spectralColors{1}};
for s = 1:2
    [~,ia,ib] = intersect(results.(whatSpecies{s}).intra.GOID,results.(whatSpecies{s}).randomReal.GOID);
    GOTableCombined = results.(whatSpecies{s}).intra(ia,:);
    GOTableCombined.FPSR_random = results.(whatSpecies{s}).randomReal.sumUnderSig(ib)/numNullSamples_surrogate;
    BF_PlotQuantiles(GOTableCombined.(theXfield),GOTableCombined.(theYfield),...
            numQuantileBins,false,false,theColors{s},false);
end
ylabel('FPSR (random null)')
xlabel(sprintf('intra-category %s',whatIntraStat))
f.Position = [1000        1117         326         221];


end
