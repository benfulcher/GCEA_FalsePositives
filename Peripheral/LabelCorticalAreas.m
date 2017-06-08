function [areaLabels,labelInd,corticalFields] = LabelCorticalAreas(regionNames)
% Give in some region names, and matches to group names for parts of the cortex
%-------------------------------------------------------------------------------

% First set up the cortical labels structure:
corticalLabels = struct();
corticalLabels.MedicalPrefrontal = {'ILA','ORBm','PL'};
corticalLabels.Lateral = {'AIv','ECT','GU','PERI','TEa','AIp','AId','VISC'};
corticalLabels.MedialAssociation = {'RSPagl','ORBvl','RSPv','ACAv','ACAd',...
                                        'RSPd','ORBl','PTLp'};
corticalLabels.MotorSomatoSensory = {'MOs','MOp','SSp-tr','SSp-ll','SSp-ul',...
                                        'SSp-m','SSp-n','SSp-bfd','SSs'};
corticalLabels.AudioVisual = {'VISam','VISl','VISal','VISpm','VISp','VISpl',...
                                        'AUDd','AUDp','AUDv','AUDpo'};
if ismember('FRP',regionNames)
    corticalLabels.Other = {'FRP'};
end

corticalLabels = structfun(@lower,corticalLabels,'UniformOutput',0);

%-------------------------------------------------------------------------------
% Then do the matching:
numRegions = length(regionNames);
regionNames = lower(regionNames);
areaLabels = cell(numRegions,1);
labelInd = zeros(numRegions,1);
corticalFields = fieldnames(corticalLabels);
for i = 1:numRegions
    whereMatch = structfun(@(x)any(ismember(regionNames{i},x)),corticalLabels);
    try
        areaLabels{i} = corticalFields{whereMatch};
    catch
        warning('Could not find %s',regionNames{i})
        keyboard
    end
    labelInd(i) = find(whereMatch);
end

end
