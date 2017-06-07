
%-------------------------------------------------------------------------------
% Connectome parameters
%-------------------------------------------------------------------------------
connectomeSource = 'Oh'; % 'Oh-cortex'
pThreshold = 0.05;
whatHemispheres = 'right';
whatWeightMeasure = 'NCD';
justCortex = true;

[A_bin,regionInfo,adjPVals] = GiveMeAdj(connectomeSource,pThreshold,true,...
                                whatWeightMeasure,whatHemispheres,justCortex);

%-------------------------------------------------------------------------------
% Gene data
%-------------------------------------------------------------------------------
[geneData,geneInfo,structInfo] = LoadMeG({'none','none'},'energy');
% Match to regions:
[~,~,isConnectome] = intersect({regionInfo.acronym},structInfo.acronym,'stable');
geneData = geneData(isConnectome,:);
structInfo = structInfo(isConnectome,:);
geneDataZ = BF_NormalizeMatrix(geneData,'zscore');

%-------------------------------------------------------------------------------
% Mean expression enrichment
%-------------------------------------------------------------------------------
gScore = nanmean(geneData,1);
% gScore = gScore(randperm(length(gScore)));
numStructs = height(structInfo);
[GOTable,geneEntrezAnnotations] = SingleEnrichment(gScore,geneInfo.entrez_id,'biological_process',[5,200],20000);

%-------------------------------------------------------------------------------
% Expression variance enrichment
%-------------------------------------------------------------------------------
gScore = nanstd(geneData,1);
% gScore = gScore(randperm(length(gScore)));
numStructs = height(structInfo);
[GOTable,geneEntrezAnnotations] = SingleEnrichment(gScore,geneInfo.entrez_id,'biological_process',[5,200],20000);

%-------------------------------------------------------------------------------
% Degree enrichment:
%-------------------------------------------------------------------------------
k = sum(A_bin,1)' + sum(A_bin,2);
gScore = zeros(height(geneInfo),1);
for i = 1:height(geneInfo)
    gScore(i) = corr(k,geneData(:,i),'type','Pearson','rows','pairwise');
end
% Plot:
[~,iy] = sort(gScore,'descend');
% [~,ix] = sort(k,'descend');
ix = 1:numStructs;
PlotGeneExpression(geneData(ix,iy),geneInfo(iy,:),structInfo(ix,:),true,k(ix))
% Enrichment:
[GOTable,geneEntrezAnnotations] = SingleEnrichment(gScore,geneInfo.entrez_id,'biological_process',[5,200],20000);

%-------------------------------------------------------------------------------
% Cerebral cortex enrichment
%-------------------------------------------------------------------------------
gScore = zeros(height(geneInfo),1);
isCTX = ismember(structInfo.divisionLabel,{'Isocortex','Olfactory Areas','Hippocampal Formation','Cortical Subplate'});
for i = 1:height(geneInfo)
    gScore(i) = ttest2(geneDataZ(isCTX,i),geneDataZ(~isCTX,i),'VarType','Unequal');
end
[GOTable,geneEntrezAnnotations] = SingleEnrichment(gScore,geneEntrezIDs,'biological_process',[5,200],20000);

%-------------------------------------------------------------------------------
% Get region-based anatomical nulls:
%-------------------------------------------------------------------------------
[GOTable,geneEntrezAnnotations] = GetFilteredGOData('biological_process',[5,200],geneInfo.entrez_id);
sizeGOCategories = cellfun(@length,geneEntrezAnnotations);
numGOCategories = height(GOTable);
numNulls = 200;
categoryScores = nan(numGOCategories,numNulls+1);
k = sum(A_bin,1)' + sum(A_bin,2);
parfor n = 1:numNulls+1
    fprintf(1,'Null %u/%u\n',n,numNulls+1);
    if n == 1
        permVector = 1:size(geneData,1);
    else
        % permVector = randperm(size(geneData,1));
        permVector = AnatomyShuffle(structInfo.divisionLabel,'twoBroad');
    end

    gScore = zeros(height(geneInfo),1);
    for i = 1:height(geneInfo)
        gScore(i) = corr(k,geneData(permVector,i),'type','Spearman','rows','pairwise');
    end

    % Record mean scores for each category:
    for j = 1:numGOCategories
        matchMe = ismember(geneInfo.entrez_id,geneEntrezAnnotations{j});
        if sum(matchMe) <= 1
            continue
        end
        categoryScores(j,n) = nanmean(gScore(matchMe));
    end
end

%-------------------------------------------------------------------------------
% Estimate p-values:
%-------------------------------------------------------------------------------
whatTail = 'right';
[meanNull,stdNull,pValsPerm,pValsZ,pValsZ_corr] = EstimatePVals(categoryScores,numNulls,whatTail);
ListCategories(geneInfo,GOTable,geneEntrezAnnotations,meanNull,pValsZ,pValsZ_corr);
