function GOTableGeneric = SurrogateEnrichmentProcess(whatSpecies,numMaps,whatSurrogate,customSurrogate)

%-------------------------------------------------------------------------------
% Check inputs:
if nargin < 1 || isempty(whatSpecies)
    whatSpecies = 'mouse';
    % whatSpecies = 'human';
end
if nargin < 2 || isempty(numMaps)
    numMaps = 10000;
end
if nargin < 3 || isempty(whatSurrogate)
    whatSurrogate = 'spatialLag';
    % whatSurrogate = 'independentSpatialShuffle';
    % whatSurrogate = 'geneShuffle';
end
if nargin < 4
    customSurrogate = '';
end

%-------------------------------------------------------------------------------
fileNameIn = sprintf('SurrogateGOTables_%u_%s_%s_%s.mat',numMaps,whatSpecies,whatSurrogate,customSurrogate);
load(fileNameIn,'GOTableGeneric','surrogatePVals');
fprintf(1,'(Data loaded from %s)\n',fileNameIn);
fprintf(1,'Enrichment of %s nulls under a %s model\n',whatSpecies,whatSurrogate);

%-------------------------------------------------------------------------------
% Compute corrected p-vals:
pValCorr = zeros(size(surrogatePVals));
for j = 1:size(surrogatePVals,2)
    pValCorr(:,j) = mafdr(surrogatePVals(:,j),'BHFDR','true');
end

%-------------------------------------------------------------------------------
% Get some statistics:
sumSig = sum(pValCorr < 0.05,2);
GOTableGeneric.sumUnderSig = sumSig;
GOTableGeneric = sortrows(GOTableGeneric,'sumUnderSig','descend');

%-------------------------------------------------------------------------------
% Estimate a p-value? NB NB NB: THIS IS WRONG!!!:
% GOTableGeneric.pValCorrNot = 1./GOTableGeneric.sumUnderSig;

%-------------------------------------------------------------------------------
display(GOTableGeneric(1:20,:))


end
