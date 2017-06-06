function permVector = AnatomyShuffle(divisionLabels)
% Shuffles regions based on anatomy

numRegions = length(divisionLabels);
groups = zeros(numRegions,1);

% CTX:
isCTX = ismember(divisionLabels,{'Isocortex','Olfactory Areas','Hippocampal Formation','Cortical Subplate'});
groups(isCTX) = 1;

% Striatum:
isStriatum = strcmp(divisionLabels,'Striatium');
groups(isStriatum) = 2;

% Pallidum Interbrain
isPallidumInterbrain = ismember(divisionLabels,{'Pallidum','Thalamus','Hypothalamus'});
groups(isPallidumInterbrain) = 3;

% Midbrain,Hindbrain
isMidHindbrain = ismember(divisionLabels,{'Midbrain','Pons','Medulla'});
groups(isMidHindbrain) = 4;

% Cerebellum:
isCerebellum = ismember(divisionLabels,{'Cerebellar Cortex','Cerebellar Nuclei'});
groups(isCerebellum) = 5;

%-------------------------------------------------------------------------------
permVector = 1:numRegions;
for i = 1:5
    theIdx = permVector(groups==i);
    permVector(groups==i) = theIdx(randperm(length(theIdx)));
end


end
