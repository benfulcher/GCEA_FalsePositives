function [categoryScores,categoryInfo] = GiveMeCategoryNullDist(whatGOID,params,numNullSamples,whatCorr,doCompute)
% Computes null distribution for a GO category
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
if nargin < 1
    whatGOID = 6099;
end
if nargin < 2
    params = GiveMeDefaultParams('mouse');
end
if nargin < 3
    numNullSamples = 100;
end
if nargin < 4
    whatCorr = 'Spearman';
end
if nargin < 5
    doCompute = false;
end
beVerbose = false;
%-------------------------------------------------------------------------------

if ~doCompute
    load('RandomNull_20000_mouse_randomMap_Spearman_mean.mat');
    
else
    %-------------------------------------------------------------------------------
    % Load in real gene-expression data for this GO category:
    [geneData,geneInfo,structInfo,categoryInfo] = GiveMeGOCategory(whatGOID,params);
    numGenesCategory = height(geneInfo);
    % PlotCategoryExpression(whatGOID,params)
    % PlotCategoryIntraCorr(whatGOID,params,whatCorr)

    %-------------------------------------------------------------------------------
    % Get random vectors from real genes to use as null spatial maps:
    params.g.humanOrMouse = sprintf('surrogate-%s',params.g.humanOrMouse);
    nullMaps = LoadMeG(params.g);

    %-------------------------------------------------------------------------------
    % Compute the distribution of gene category scores for correlation with the null maps:
    categoryScores = nan(numNullSamples,1);
    for i = 1:numNullSamples
        if beVerbose
            fprintf(1,'%u/%u\n',i,numNullSamples);
        end
        geneScores = zeros(numGenesCategory,1);
        for j = 1:numGenesCategory
            geneScores(j) = corr(nullMaps(:,i),geneData(:,j),'type',whatCorr,'rows','pairwise');
        end
        categoryScores(i) = nanmean(geneScores);
    end
end

end
