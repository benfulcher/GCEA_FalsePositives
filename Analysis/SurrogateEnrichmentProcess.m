
% whatSpecies = 'mouse';
whatSpecies = 'human';


switch whatSpecies
case 'mouse'
    load('SurrogateGOTables_1000_mouse_spatialLag.mat','GOTableGeneric','surrogatePVals');
case 'human'
    load('SurrogateGOTables_1000_human_spatialLag.mat','GOTableGeneric','surrogatePVals')
end


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
