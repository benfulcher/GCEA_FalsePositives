%-------------------------------------------------------------------------------
% Annotate intra-category correlations for different gene-expression datasets
%-------------------------------------------------------------------------------
% Parameters:
whatSpecies = 'mouse';
whatSurrogate = 'geneShuffle';
numSamples = 5000;

%-------------------------------------------------------------------------------
% Set default parameters:
params = GiveMeDefaultParams(whatSpecies);
params.g.whatSurrogate = whatSurrogate;

%-------------------------------------------------------------------------------
% Load in expression data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);

%-------------------------------------------------------------------------------
% Compute intra-category correlations:
resultsTable = AnnotateIntraCorrelations(params,[],whatSpecies);
numGOCategories = height(resultsTable);

%-------------------------------------------------------------------------------
% Compute a null distribution for each category size:
uniqueSizes = unique(resultsTable.size);
numSizes = length(uniqueSizes);
numGenes = height(geneInfo);

nullDistributionRaw = zeros(numSizes,numSamples);
nullDistributionAbs = zeros(numSizes,numSamples);
switch whatSurrogate
case 'spatialShuffle'
    % Spatial shuffle case:
    fprintf(1,'(INDEPENDENT) SPATIAL SHUFFLE!!\n');
    params.g.humanOrMouse = 'surrogate-mouse';
    params.g.whatSurrogate = 'spatialShuffle';
    parfor j = 1:numSamples
        fprintf(1,'Sample %u/%u\n',j,numSamples);
        % Get indpendently-shuffled data:
        [geneData,geneInfo,structInfo] = LoadMeG(params.g);
        for i = 1:numSizes
            shuffledData = geneData(:,1:uniqueSizes(i));
            [nullDistributionRaw(i,j),nullDistributionAbs(i,j)] = IntraCorrelationScore(shuffledData);
        end
    end
case 'geneShuffle'
    % Shuffle genes randomly:
    fprintf(1,'GENE SHUFFLE!!\n');
    nullDistributionRaw = zeros(numSizes,numSamples);
    nullDistributionAbs = zeros(numSizes,numSamples);
    parfor j = 1:numSamples
        fprintf(1,'Sample %u/%u\n',j,numSamples);
        % Get indpendently shuffled data:
        geneDataShuffle = geneData(:,randperm(numGenes));
        for i = 1:numSizes
            shuffledData = geneDataShuffle(:,1:uniqueSizes(i));
            [nullDistributionRaw(i,j),nullDistributionAbs(i,j)] = IntraCorrelationScore(shuffledData);
        end
    end
end

%-------------------------------------------------------------------------------
theField = 'mouse_abs';
theNullDistribution = nullDistributionAbs;

%-------------------------------------------------------------------------------
% Plot:
f = figure('color','w'); hold('on')
h1 = histogram(resultsTable.(theField),'normalization','pdf');
h2 = histogram(theNullDistribution(1,:),'normalization','pdf');
h3 = histogram(theNullDistribution(end,:),'normalization','pdf');
legend([h1,h2,h3],{'realData',sprintf('nullSize%u',uniqueSizes(1)),sprintf('nullSize%u',uniqueSizes(end))});

%-------------------------------------------------------------------------------
% Compute p-values (bigger scores are better)
pVal = nan(numGOCategories,1);
parfor i = 1:numGOCategories
    categoryScore = resultsTable.(theField)(i);
    if ~isnan(categoryScore)
        nullForSize = theNullDistribution(resultsTable.size(i)==uniqueSizes,:);
        pVal(i) = mean(categoryScore < nullForSize);
    end
end

%-------------------------------------------------------------------------------
% FDR correct:
pValCorr = mafdr(pVal,'BHFDR','true');

%-------------------------------------------------------------------------------
% Update the GO table:
resultsTable.pVal = pVal;
resultsTable.pValCorr = pValCorr;

%-------------------------------------------------------------------------------
% Sort:
resultsTable = sortrows(resultsTable,{'pVal','pValCorr'},{'ascend','ascend'});
