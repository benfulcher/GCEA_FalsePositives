function GOTableGeneric = SurrogateEnrichmentProcess(params,doDisplay)
% Process precomputed FPSR results for analysis
%-------------------------------------------------------------------------------
% Check inputs:
if nargin < 1
    params = GiveMeDefaultParams('mouse');
end
if nargin < 2
    doDisplay = true;
end

%-------------------------------------------------------------------------------
fileNameFPSR = GiveMeFPSRFileName(params);
load(fileNameFPSR,'GOTableGeneric')
fprintf(1,'(Data loaded from %s)\n',fileNameFPSR);
if params.nulls.permTestP
    load(fileNameFPSR,'surrogatePValsPerm');
    surrogatePVals = surrogatePValsPerm;
    clear('surrogatePValsPerm');
    fprintf(1,'Using permutation test p-values\n');
else
    load(fileNameFPSR,'surrogatePValsZ');
    surrogatePVals = surrogatePValsZ;
    clear('surrogatePValsZ');
    fprintf(1,'Using Gaussian-approx permutation test p-values\n');
end
fprintf(1,'Enrichment of %s nulls under a %s model\n',...
                        params.humanOrMouse,params.g.whatSurrogate);

%-------------------------------------------------------------------------------
% Compute corrected p-vals:
pValCorr = zeros(size(surrogatePVals));
for j = 1:size(surrogatePVals,2)
    pValCorr(:,j) = mafdr(surrogatePVals(:,j),'BHFDR','true');
end

%-------------------------------------------------------------------------------
% Get some statistics:
sumSig = sum(pValCorr < params.e.sigThresh,2);
GOTableGeneric.sumUnderSig = sumSig;
GOTableGeneric = sortrows(GOTableGeneric,'sumUnderSig','descend');

%-------------------------------------------------------------------------------
% Show the top XXX cateogries:
if doDisplay
    display(GOTableGeneric(1:20,:))
end

end
