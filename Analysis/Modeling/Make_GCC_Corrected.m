function Make_GCC_Corrected(removeDeadEnd)
%-------------------------------------------------------------------------------
% Quick function to try to extract a smaller file containing corrected GCC scores
% (which can maybe be run locally)
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Set Parameters:
whatGeneData = 'norm'; % 'norm', 'raw', 'zscore'
EnergyOrDensity = 'energy'; % 'energy', 'density'
theCorrection = 'distance'; % 'intraInter','divII','distance','none','cortex','intralinear'
theFitType = 'expFitAll'; %,'exp_1_0', 'expFitAll', 'linear','decayEta','exp_1_0','expFitLinked'
baselineLinks = 'all'; % 'all', 'connected' (do the Ghat correction for what set of links?)

%-------------------------------------------------------------------------------
% Load the GCC scores from file
extraText = sprintf('_%s_%s',whatGeneData,EnergyOrDensity);
fprintf(1,'Loading in ggBlocks and Ghats data for this %s/%s analysis... ',...
                            whatGeneData,EnergyOrDensity);
load(['ggBlocks',extraText,'.mat'],'ggBlocks_raw','theGeneStruct'); % this is 5.7GB :-/
numGenes = size(ggBlocks_raw,1);
load(['Ghats',extraText,'.mat'],'Ghat')
fprintf(1,' Done.\n');

%-------------------------------------------------------------------------------
% Apply a global correction to all GCC scores:
fprintf(1,'Correcting raw ggBlocks on the basis of %s/%s/%s...',...
                            theCorrection,theFitType,baselineLinks);
ggBlockscorr = zeros(size(ggBlocks_raw));
Gcorrection = Ghat.(theCorrection).(theFitType).(baselineLinks);
for i = 1:numGenes
    rawData = squeeze(ggBlocks_raw(i,:,:));
    ggBlockscorr(i,:,:) = rawData - Gcorrection;
end
fprintf(1,'Done!!!\n');
clear('Ghat','ggBlocks_raw');

%-------------------------------------------------------------------------------
% Get the elements where links of binary A at 0.05 exist:
C = load('Mouse_Connectivity_Data.mat'); % C stands for connectome data
% Get the adjacency matrix and process it:
A = logical(GiveMeAdj(C,'binary','ipsi',0,0.05));

if removeDeadEnd
    k_out = sum(A,2); % out degree
    deadEndNodes = (k_out==0);
    A = A(~deadEndNodes,~deadEndNodes);
    ggBlockscorr = ggBlockscorr(:,~deadEndNodes,~deadEndNodes);
    fprintf(1,'Removed %u dead end nodes\n',sum(deadEndNodes));
end

%-------------------------------------------------------------------------------
numLinks = sum(A(:));
ggCorr_A_p005 = zeros(numGenes,numLinks);
for i = 1:numGenes
    GCC_i = squeeze(ggBlockscorr(i,:,:));
    ggCorr_A_p005(i,:) = GCC_i(A);
end


if removeDeadEnd
    fileName = 'GCC_A_p005_deadEnd.mat';
else
    fileName = 'GCC_A_p005.mat';
end

save(fileName,'ggCorr_A_p005','theGeneStruct','whatGeneData','EnergyOrDensity',...
                'theCorrection','theFitType','baselineLinks','removeDeadEnd');

end
