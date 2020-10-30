% Computation for all case studies -> AllCaseStudies.mat

%===============================================================================
% DEGREE
%===============================================================================
nodeMetrics = {'degree','betweenness'};
numNodeMetrics = length(nodeMetrics);

% Precomputing conventional results:
for n = 1:numNodeMetrics
    params = GiveMeDefaultParams('mouse','cortex');
    NodeSimpleEnrichment(params,nodeMetrics{n},true);

    params = GiveMeDefaultParams('mouse','all');
    NodeSimpleEnrichment(params,nodeMetrics{n},true);

    params = GiveMeDefaultParams('human','cortex');
    NodeSimpleEnrichment(params,nodeMetrics{n},true);
end

%-------------------------------------------------------------------------------
% Degree/Betweenness:
for n = 1:numNodeMetrics
    % Mouse brain
    T = CaseStudyResults('mouse','all',nodeMetrics{n});
    T = sortrows(T,'pValZRandomGene');
    resultsTables.(sprintf('mouseBrain%s',nodeMetrics{n})) = T;
    T(1:20,:)

    % Mouse cortex
    T = CaseStudyResults('mouse','cortex',nodeMetrics{n});
    T = sortrows(T,'pValZRandomGene');
    resultsTables.(sprintf('mouseCortex%s',nodeMetrics{n})) = T;
    T(1:20,:)

    % Human cortex
    T = CaseStudyResults('human','cortex',nodeMetrics{n});
    T = sortrows(T,'pValZRandomGene');
    resultsTables.(sprintf('humanCortex%s',nodeMetrics{n})) = T;
    T(1:20,:)
end

%===============================================================================
% CELL DENSITIES: mouse cortex/brain
%===============================================================================
cellTypes = {'excitatory','inhibitory','oligodendrocytes','glia','astrocytes',...
                        'microglia','neurons','PV','SST','VIP'};
numCellTypes = length(cellTypes);

%-------------------------------------------------------------------------------
% Precomputing conventional results (mouse cortex/brain):
params_cortex = GiveMeDefaultParams('mouse','cortex');
params_brain = GiveMeDefaultParams('mouse','all');
for c = 1:numCellTypes
    fprintf(1,'Precomputing random-gene results for %s\n',cellTypes{c});
    % NodeSimpleEnrichment(params_cortex,cellTypes{c},true);
    NodeSimpleEnrichment(params_brain,cellTypes{c},true);
end

%-------------------------------------------------------------------------------
% Now compare the different null models:
for c = 1:numCellTypes
    % ---Mouse cortex---
    % T = CaseStudyResults('mouse','cortex',cellTypes{i});
    % T = sortrows(T,'pValZRandomGene');
    % T(1:20,:)
    % resultsTables.(sprintf('mouseCortex_%s',cellTypes{i})) = T;
    % ---Mouse brain---
    T = CaseStudyResults('mouse','all',cellTypes{c});
    T = sortrows(T,'pValZRandomGene');
    T(1:20,:)
    resultsTables.(sprintf('mouseBrain_%s',cellTypes{c})) = T;
end

%-------------------------------------------------------------------------------
% Save to .mat file to avoid need for further computation
fprintf(1,'Saving all results to file\n');
save('AllCaseStudies.mat','resultsTables')

%-------------------------------------------------------------------------------
% List out top hits:
ProcessCaseStudy(resultsTables)
