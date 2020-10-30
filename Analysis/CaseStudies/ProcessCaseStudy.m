function [resultsTables,numSig] = ProcessCaseStudy(resultsTables,justCortex,sigThresh,displayToScreen)
% Output a summary of all results to text file:

%-------------------------------------------------------------------------------
if nargin < 2
    justCortex = true;
end
if nargin < 3
    sigThresh = 0.05;
end
if nargin < 4
    displayToScreen = true;
end

%-------------------------------------------------------------------------------
% Filter out whole-mouse-brain results?:
if justCortex
    theAnalyses = fieldnames(resultsTables);
    isMouseBrain = ~cellfun(@isempty,regexp(theAnalyses,'mouseBrain'));
    resultsTables = rmfield(resultsTables,theAnalyses(isMouseBrain));
end
theAnalyses = fieldnames(resultsTables);
numAnalyses = length(theAnalyses);

%-------------------------------------------------------------------------------
% Display to screen
maxShow = 10;
nullNames = {'Random Gene','SBP-Random','SBP-Spatial'};
nullsRaw = {'pValZRandomGene','pValZRandomMap','pValZSpatialLag'};
nullsCorr = {'pValZCorrRandomGene','pValZCorrRandomMap','pValZCorrSpatialLag'};
numNulls = length(nullsRaw);
numSig = nan(numAnalyses,numNulls);
for i = 1:numAnalyses
    if displayToScreen
        fprintf(1,'\n\n----------%s------------\n\n',theAnalyses{i});
    end

    for k = 1:numNulls
        isSignificant = (resultsTables.(theAnalyses{i}).(nullsCorr{k}) < sigThresh);
        numSig(i,k) = sum(isSignificant);
        if displayToScreen
            fprintf(1,'%s null: %u significant (Gaussian-approx).\n\n',nullNames{k},numSig);
            [~,ix] = sort(resultsTables.(theAnalyses{i}).(nullsRaw{k}),'ascend');
            for j = 1:min(numSig(i,k),maxShow)
                fprintf(1,'%s (%g)\n',resultsTables.(theAnalyses{i}).GOName{ix(j)},...
                                resultsTables.(theAnalyses{i}).(nullsCorr{k})(ix(j)));
            end
            fprintf(1,'\n');
        end
    end
end

end
