function fileNameFPSR = GiveMeFPSRFileName(params)
% Get filename for a given FPSR analysis
%-------------------------------------------------------------------------------

fileNameFPSR = sprintf('SurrogateGOTables_%u_%s_%s_%s.mat',...
                params.nulls.numNullsCFPR,...
                params.humanOrMouse,...
                params.g.whatSurrogate,...
                params.nulls.customShuffle);

fileNameFPSR = fullfile('SurrogateEnrichment',fileNameFPSR);

end