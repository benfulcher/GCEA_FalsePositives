
% whatSpecies = 'mouse';
whatSpecies = 'human';

% whatSurrogate = 'spatialLag';
whatSurrogate = 'spatialShuffle';

%-------------------------------------------------------------------------------

theMatFile = sprintf('SurrogateGOTables_1000_%s_%s.mat',whatSpecies,whatSurrogate);
load(theMatFile,'GOTableGeneric','surrogatePVals');

fprintf(1,'Enrichment of %s nulls under a %s model\n',whatSpecies,whatSurrogate);

%-------------------------------------------------------------------------------
% Compute corrected p-vals:
pValCorr = zeros(size(surrogatePVals));
for j = 1:size(surrogatePVals,2)
    pValCorr(:,j) = mafdr(surrogatePVals(:,j),'BHFDR','true');
end

%-------------------------------------------------------------------------------
% Get some statistics:
% meanPCorr = mean(pValCorr,2);
sumUnderSig = sum(pValCorr < 0.05,2);
GOTableGeneric.sumUnderSig = sumUnderSig;

GOTableGeneric = sortrows(GOTableGeneric,'sumUnderSig','descend');

%-------------------------------------------------------------------------------
display(GOTableGeneric(1:20,:))
