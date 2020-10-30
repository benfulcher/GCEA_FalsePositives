% Load results from running AllCaseStudies:
fprintf(1,'Loading case study results to file\n');
load('AllCaseStudies.mat','resultsTables')

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
