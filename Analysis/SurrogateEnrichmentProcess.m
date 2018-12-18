function GOTableGeneric = SurrogateEnrichmentProcess(whatSpecies,whatSurrogate)
%
%-------------------------------------------------------------------------------

% Check inputs:
if nargin < 1 || isempty(whatSpecies)
    whatSpecies = 'mouse';
    % whatSpecies = 'human';
end
if nargin < 2 || isempty(whatSurrogate)
    whatSurrogate = 'spatialLag';
    % whatSurrogate = 'spatialShuffle';
    % whatSurrogate = 'geneShuffle';
end

%-------------------------------------------------------------------------------
theMatFile = sprintf('SurrogateGOTables_1000_%s_%s.mat',whatSpecies,whatSurrogate);
load(theMatFile,'GOTableGeneric','surrogatePVals');
fprintf(1,'(Data loaded from %s)\n',theMatFile);
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
GOTableGeneric.pValCorr = 1./GOTableGeneric.sumUnderSig;

%-------------------------------------------------------------------------------
% display(GOTableGeneric(1:20,:))

end
