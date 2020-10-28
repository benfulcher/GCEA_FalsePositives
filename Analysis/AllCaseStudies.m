
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
% Output a summary of all results to text file:
theAnalyses = fieldnames(resultsTables);
numAnalyses = length(theAnalyses);
sigThresh = 0.05;
maxShow = 10;
nullNames = {'Random Gene','SBP-Random','SBP-Spatial'};
nullsRaw = {'pValZRandomGene','pValZRandomMap','pValZSpatialLag'};
nullsCorr = {'pValZCorrRandomGene','pValZCorrRandomMap','pValZCorrSpatialLag'};
numNulls = length(nullsRaw);
numSig = nan(numAnalyses,numNulls);
for i = 1:numAnalyses
    fprintf(1,'\n\n----------%s------------\n\n',theAnalyses{i});

    for k = 1:numNulls
        isSignificant = (resultsTables.(theAnalyses{i}).(nullsCorr{k}) < sigThresh);
        numSig(i,k) = sum(isSignificant);
        fprintf(1,'%s null: %u significant (Gaussian-approx).\n\n',nullNames{k},numSig);
        [~,ix] = sort(resultsTables.(theAnalyses{i}).(nullsRaw{k}),'ascend');
        for j = 1:min(numSig(i,k),maxShow)
            fprintf(1,'%s (%g)\n',resultsTables.(theAnalyses{i}).GOName{ix(j)},...
                            resultsTables.(theAnalyses{i}).(nullsCorr{k})(ix(j)));
        end
        fprintf(1,'\n');
    end
end

%-------------------------------------------------------------------------------
% Extras:
%-------------------------------------------------------------------------------
% Check whether sensible inhibitory GO categories are coming up strongly with inhibitory
% cell density...?

% Search for GABA?
isGABA = find(cellfun(@(x)~isempty(x),regexp(resultsTables.mouseCortex_inhibitory.GOName,'GABA')));
display(resultsTables.mouseCortex_inhibitory(isGABA,:))

% Check the full histogram
f = figure('color','w');
histogram(-log10(resultsTables.mouseCortex_inhibitory.pValZRandomGene),'numBins',50)
hold('on');

% Check the actual correlation:
params = GiveMeDefaultParams('mouse','cortex');
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
[phenotypeVector,ia,regionNamesOrdered] = MatchMeCellDensity(structInfo,'inhibitory');
structInfo = structInfo(ia,:);
geneData = geneData(ia,:);

GOTable = GiveMeGOData(params,geneInfo.entrez_id);
isGABA = find(cellfun(@(x)~isempty(x),regexp(GOTable.GOName,'GABA')));

theGenesEntrez = GOTable.annotations{isGABA(2)};
[entrezMatched,ia,ib] = intersect(theGenesEntrez,geneInfo.entrez_id);
matchMe = find(ismember(geneInfo.entrez_id,entrezMatched));

[~,ib] = sort(phenotypeVector,'descend');

f = figure('color','w');
ax = subplot(1,6,1:5)
geneDataCategory = geneData(ib,matchMe);
geneInfoCategory = geneInfo(matchMe,:);
ord_col = BF_ClusterReorder(geneDataCategory','corr','average');
BF_imagesc(BF_NormalizeMatrix(geneDataCategory(:,ord_col),'scaledRobustSigmoid'))
ax.XTick = 1:length(matchMe);
ax.XTickLabel = geneInfoCategory.acronym(ord_col);
ax.XTickLabelRotation = 90;
title(GOTable.GOName{isGABA(2)})
ax.YTick = 1:height(structInfo);
ax.YTickLabel = structInfo.acronym(ib);
% imagesc(geneData(:,matchMe))
ax = subplot(1,6,6);
imagesc(BF_NormalizeMatrix(phenotypeVector(ib),'scaledRobustSigmoid'))
ax.XTick = 1;
ax.XTickLabel = {'inhibitory density'};
ax.YTick = [];
colormap([flipud(BF_getcmap('blues',9));1,1,1;BF_getcmap('reds',9)])

%-------------------------------------------------------------------------------
% Are correlations with inhibitory density generally higher across genes than to
% other phenotypes?
%-------------------------------------------------------------------------------
enrichWhat = {'excitatory','inhibitory','neurons'};
corrType = 'Spearman';
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
numGenes = height(geneInfo);
gScores = cell(length(enrichWhat),1);
for e = 1:length(enrichWhat)
    fprintf(1,'%s\n',enrichWhat{e});
    [phenotypeVector,ia] = MatchMeCellDensity(structInfo,enrichWhat{e});
    geneData_e = geneData(ia,:);
    gScore = zeros(numGenes,1);
    for g = 1:numGenes
        gScores{e}(g) = corr(phenotypeVector,geneData_e(:,g),'type',corrType,'rows','pairwise');
    end
end

f = figure('color','w');
hold('on')
for e = 1:3
    histogram(gScores{e},'Normalization','CountDensity');
end
legend(enrichWhat)
xlabel('Correlation')
ylabel('Count density')

%===============================================================================
% Now, how can we nicely show the results?
%===============================================================================
justCortex = true;
showZero = false;
if justCortex
    isMouseBrain = ~cellfun(@isempty,regexp(theAnalyses,'mouseBrain'));
    theAnalysesShow = theAnalyses(~isMouseBrain);
    numSigShow = numSig(~isMouseBrain,:);
else
    theAnalysesShow = theAnalyses;
    numSigShow = numSig;
end
if ~showZero
    myFilter = sum(numSigShow,2)>0;
    numSigShow = numSigShow(myFilter,:);
    theAnalysesShow = theAnalysesShow(myFilter);
end
f = figure('color','w');
ax = gca();
bar(numSigShow)
ax.XTick = 1:size(numSigShow,1);
ax.XTickLabel = theAnalysesShow;
ax.XTickLabelRotation = 30;
ax.TickLabelInterpreter = 'None';
legend('random-Gene','SBP-rand','SBP-spatial')
title('Significant enrichment across phenotypes')

