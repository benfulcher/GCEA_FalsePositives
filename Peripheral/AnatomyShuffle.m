function permVector = AnatomyShuffle(divisionLabels,whatLevel,beVocal)
% Shuffles regions based on anatomy

if nargin < 2
    whatLevel = 'fiveByEye';
end
if nargin < 3
    beVocal = true;
end

numRegions = length(divisionLabels);
groups = zeros(numRegions,1);

switch whatLevel
case 'fiveByEye'
    if beVocal
        fprintf(1,'SHUFFLING USING FIVE ANATOMICAL GROUPS\n');
    end
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
case 'twoIsocortex'
    if beVocal
        fprintf(1,'SHUFFLING SEPARATELY WITHIN ISOCORTEX & NON-ISOCORTEX\n');
    end
    isCTX = strcmp(divisionLabels,'Isocortex');
    groups(isCTX) = 1;
    groups(~isCTX) = 2;
case 'twoBroad'
    if beVocal
        fprintf(1,'SHUFFLING USING TWO ANATOMICAL GROUPS\n');
    end
    isCTX = ismember(divisionLabels,{'Isocortex','Olfactory Areas','Hippocampal Formation','Cortical Subplate'});
    groups(isCTX) = 1;
    groups(~isCTX) = 2;
end

%-------------------------------------------------------------------------------
permVector = 1:numRegions;
for i = 1:max(groups)
    theIdx = permVector(groups==i);
    permVector(groups==i) = theIdx(randperm(length(theIdx)));
end


end
