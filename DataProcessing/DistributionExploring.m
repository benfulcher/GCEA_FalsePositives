% Idea is to take a look at the distribution of the data
% (and examine the appropriateness of normalization)
%-------------------------------------------------------------------------------

doNormalize = false;
energyOrDensity = 'energy';
normalizationSettings = {'log10','zscore'}; % {geneNormalization,regionNormalization}

%-------------------------------------------------------------------------------
[GeneStruct,GeneExpData] = LoadMeG(normalizationSettings,energyOrDensity);
[numRegions,numGenes] = size(GeneExpData);

%-------------------------------------------------------------------------------
% Normalize:
% switch whatNormalizationGenes
% case 'log'
%     GeneExpDataN = GeneExpData;
%     GeneExpDataN(GeneExpDataN==0) = NaN;
%     GeneExpDataN = log10(GeneExpData);
% otherwise
%     GeneExpDataN = BF_NormalizeMatrix(GeneExpData,whatNormalizationGenes);
% end
% GeneExpDataN = BF_NormalizeMatrix(GeneExpDataN',whatNormalizationRegions)';

%-------------------------------------------------------------------------------
% Distribution of genes across brain regions for a sample of genes:
sampleSize = 25;
plotHow = 'scatterhist';
mySample = randi(numGenes,[sampleSize,1]);
f = figure('color','w');
for i = 1:sampleSize
    subplot(5,5,i); hold on
    switch plotHow
    case 'scatterhist'
        h = scatterhist(GeneExpData(:,mySample(i)),GeneExpDataN(:,mySample(i)));
        h(1).XLabel.String = 'raw';
        h(1).YLabel.String = sprintf('%s-%s',whatNormalizationGenes,whatNormalizationRegions);
        input('helloBois')
    case 'scatter'
        plot(GeneExpData(:,mySample(i)),GeneExpDataN(:,mySample(i)))
        xlabel('raw')
        ylabel(sprintf('%s-%s',whatNormalizationGenes,whatNormalizationRegions))
    case 'histograms'
        histogram(BF_NormalizeMatrix(GeneExpData(:,mySample(i)),'maxmin'));
        histogram(BF_NormalizeMatrix(GeneExpDataN(:,mySample(i)),'maxmin'));
        if i==1
            legend('raw',whatNormalization)
        end
    end
    title(GeneStruct(mySample(i)).gene_acronym)
end

%-------------------------------------------------------------------------------
% Distribution of genes in a given brain region:
sampleSize = 25;
mySample = randi(numRegions,[sampleSize,1]);
f = figure('color','w');
for i = 1:sampleSize
    subplot(5,5,i); hold on
    g = GeneExpData(mySample(i),:)';
    histogram(BF_NormalizeMatrix(g,'maxmin'));
    switch whatNormalization
    case 'log'
        normData = log10(g);
    otherwise
        normData = BF_NormalizeMatrix(g,whatNormalization);
    end
    histogram(BF_NormalizeMatrix(normData,'maxmin'));
    title(GeneStruct(i).gene_acronym)
    if i==1
        legend('raw',whatNormalization)
    end
end

%-------------------------------------------------------------------------------
% Distribution of expression across brain regions:
