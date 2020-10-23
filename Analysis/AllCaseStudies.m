
%===============================================================================
% DEGREE
%===============================================================================

% Precomputing conventional results:
params = GiveMeDefaultParams('mouse','cortex');
[GOTable,gScore] = NodeSimpleEnrichment(params,'degree',true);

params = GiveMeDefaultParams('mouse','all');
[GOTable,gScore] = NodeSimpleEnrichment(params,'degree',true);

params = GiveMeDefaultParams('human','cortex');
[GOTable,gScore] = NodeSimpleEnrichment(params,'degree',true);

%-------------------------------------------------------------------------------
% Degree: whole mouse brain
T = CaseStudyResults('mouse','all','degree');
T = sortrows(T,'pValZRandomGene');
resultsTables.mouseBrainDegree = T;
T(1:20,:)

%-------------------------------------------------------------------------------
% Degree: mouse cortex
T = CaseStudyResults('mouse','cortex','degree');
T = sortrows(T,'pValZRandomGene');
resultsTables.mouseCortexDegree = T;
T(1:20,:)

%-------------------------------------------------------------------------------
% Degree: human cortex brain
T = CaseStudyResults('human','cortex','degree');
T = sortrows(T,'pValZRandomGene');
resultsTables.humanCortexDegree = T;
T(1:20,:)

%===============================================================================
% CELL DENSITIES: mouse cortex
%===============================================================================

%-------------------------------------------------------------------------------
% Precomputing conventional results:
params = GiveMeDefaultParams('mouse','cortex');
cellTypes = {'excitatory','inhibitory','oligodendrocytes','glia','astrocytes',...
                        'microglia','neurons','PV','SST','VIP'};
numCellTypes = length(cellTypes);
for i = 1:numCellTypes
    fprintf(1,'Precomputing random-gene results for %s\n',cellTypes{i});
    [GOTable,gScore] = NodeSimpleEnrichment(params,cellTypes{i},true);
end
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Now compare the different null models:
for i = 1:numCellTypes
    T = CaseStudyResults('mouse','cortex',cellTypes{i});
    T = sortrows(T,'pValZRandomGene');
    T(1:20,:)
    resultsTables.(sprintf('mouseCortex_%s',cellTypes{i})) = T;
end

%-------------------------------------------------------------------------------
fprintf(1,'Saving all results to file\n');
save('AllCaseStudies.mat','resultsTables')

%-------------------------------------------------------------------------------
theAnalyses = fieldnames(resultsTables);
numAnalyses = length(theAnalyses);
sigThresh = 0.05;
maxShow = 10;
nullNames = {'Random Gene','SBP-Random','SBP-Spatial'};
nullsRaw = {'pValZRandomGene','pValZRandomMap','pValZSpatialLag'};
nullsCorr = {'pValZCorrRandomGene','pValZCorrRandomMap','pValZCorrSpatialLag'};
numNulls = length(nullsRaw);
for i = 1:numAnalyses
    fprintf(1,'\n\n----------%s------------\n\n',theAnalyses{i});

    for k = 1:numNulls
        isSignificant = (resultsTables.(theAnalyses{i}).(nullsCorr{k}) < sigThresh);
        numSig = sum(isSignificant);
        fprintf(1,'%s null: %u significant (Gaussian-approx).\n\n',nullNames{k},numSig);
        [~,ix] = sort(resultsTables.(theAnalyses{i}).(nullsRaw{k}),'ascend');
        for j = 1:min(numSig,maxShow)
            fprintf(1,'%s (%g)\n',resultsTables.(theAnalyses{i}).GOName{ix(j)},...
                            resultsTables.(theAnalyses{i}).(nullsCorr{k})(ix(j)));
        end
        fprintf(1,'\n');
    end
end
