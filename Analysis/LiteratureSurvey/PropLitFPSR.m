function PropLitFPSR(whatSpecies,whatStruct,makeEquiprobable,makeNewFigure)
% Proportion of literature results that sit in each FPSR bin
%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    whatSpecies = 'human';
end
if nargin < 2
    whatStruct = 'cortex';
end
if nargin < 3
    makeEquiprobable = false;
end
if nargin < 4
    makeNewFigure = true;
end
%-------------------------------------------------------------------------------

whatSurrogate = {'randomMap','spatialLag'};
params = cell(2,1);
LitTable = cell(2,1);
GOTable_FPSR = cell(2,1);
for s = 1:2
    params{s} = GiveMeDefaultParams(whatSpecies,whatStruct);
    params{s}.g.whatSurrogate = whatSurrogate{s};

    %-------------------------------------------------------------------------------
    % Retrieve information about how literature results are distributed across categories:
    LitTable{s} = MakeLiteratureTable(whatSpecies,params{s}.e.sigThresh);

    %-------------------------------------------------------------------------------
    % Now we'll get the FPSR data:
    GOTable_FPSR{s} = SurrogateEnrichmentProcess(params{s},false);
end

%-------------------------------------------------------------------------------
% PLOT:
% Proportion of literature appearing in each bin
if makeNewFigure
    f = figure('color','w');
end
% hold('on')
numBins = 8;
doLinear = false;
theColors = GiveMeColors('nullModels');
for s = 1:2
    myGOTable = GOTable_FPSR{s};
    % ---Filter out zero-FPSR categories---
    myGOTable = myGOTable(myGOTable.sumUnderSig > 0,:);
    if makeEquiprobable
        binEdges = arrayfun(@(x)quantile(myGOTable.sumUnderSig,x),linspace(0,1,numBins+1));
        binEdges(1) = 1; % exclude categories with zero FPSR
        binEdges(end) = binEdges(end) + eps; % make sure all data included in final bin
        binEdges = log10(binEdges); % not sure why I'm doing this
    else
        binEdges = linspace(0,max(log10(myGOTable.sumUnderSig))+eps,numBins+1);
    end
    binMeans = mean([binEdges(1:end-1);binEdges(2:end)],1);

    hasBeenReportedInBin = zeros(numBins,1);
    TableInBin = @(x,Gtable) Gtable(log10(Gtable.sumUnderSig) >= binEdges(x) & log10(Gtable.sumUnderSig) < binEdges(x+1),:);
    for i = 1:numBins
        tableInBini = TableInBin(i,GOTable_FPSR{s});
        hasBeenReportedInBin(i) = mean(ismember(LitTable{s}.GOID,tableInBini.GOID));
    end
    switch whatSurrogate{s}
    case 'randomMap'
        theColor = theColors{2};
    case 'spatialLag'
        theColor = theColors{3};
    end
    if doLinear
        plot(10.^binMeans/10000,hasBeenReportedInBin,'o-','color',theColor);
    else
        semilogx(10.^binMeans/10000,hasBeenReportedInBin,'o-','color',theColor);
    end
    hold('on')
end
xlabel('Category FPR')
ylabel({'Proportion of all literature-';'reported categories in bin'})
legend(whatSurrogate,'Location','NorthWest')
title(whatSpecies)

end
