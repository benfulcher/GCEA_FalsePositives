
params = GiveMeDefaultParams('mouse');

%===============================================================================
% CORTICAL EXPRESSION
%===============================================================================
% Gene data:
[geneData,geneInfo,structInfo_gene] = LoadMeG(params.g);
% numStructs = height(structInfo);
numGenes = height(geneInfo);

% Degree data:
[k,structInfo_k] = ComputeDegree(params,true);

% Check for structure match:
[~,ia,ib] = intersect(structInfo_gene.acronym,structInfo_k.acronym,'stable');
assert(all(diff(ia)==1) & all(diff(ia)==1))
structInfo = structInfo_gene(ia,:);

isCTX = ismember(structInfo.divisionLabel,'Isocortex');

% fprintf(1,'~~RANKSUM P-VALUE GENE SCORES~~\n');
gScore = zeros(numGenes,2);
for i = 1:numGenes
    % cortex/non-cortex expression
    gScore(i,1) = ranksum(geneData(isCTX,i),geneData(~isCTX,i));
    % degree correlation
    gScore(i,2) = corr(geneData(:,i),k,'rows','pairwise');
end

%-------------------------------------------------------------------------------
% Plot
f = figure('color','w');
xData = -log10(gScore(:,1));
yData = gScore(:,2);
plot(xData,yData,'.k')
r = corr(xData,yData,'type','Spearman','rows','pairwise');
title(sprintf('%u genes, r = %.3f',numGenes,r))
axis('square')
xlabel('-log10(p) [isocortex]')
ylabel('corr(k)')
f.Position = [1000,1027,432,311];


% % How many are individually significant?:
% pCorr = mafdr(gScore,'BHFDR','true');
% fprintf(1,'%u/%u genes have corrected p < 0.05\n',sum(pCorr < 0.05),length(gScore));
%
% % Transform p-values to scores (bigger is better)
% fprintf(1,'-log10 p-values -> gene scores\n');
% gScore = -log10(gScore);
%===============================================================================

%===============================================================================
% CORTICAL EXPRESSION
%===============================================================================
