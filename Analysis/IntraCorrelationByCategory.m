function resultsTable = IntraCorrelationByCategory(params,whatSurrogate,numSamples,pValsFromWhat)
% Annotate intra-category coexpression for different gene-expression datasets
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check input parameters:
if nargin < 1
    whatSpecies = 'mouse';
    params = GiveMeDefaultParams(whatSpecies);
end
if nargin < 2
    whatSurrogate = 'geneShuffle';
end
if nargin < 3
    numSamples = 500;
end
if nargin < 4
    pValsFromWhat = 'VE1';
end

%-------------------------------------------------------------------------------
% Set default parameters:
params.g.whatSurrogate = whatSurrogate;

%-------------------------------------------------------------------------------
% Load in expression data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);

%-------------------------------------------------------------------------------
% Compute intra-category correlations:
resultsTable = AnnotateIntraCorrelations(params,[]);
numGOCategories = height(resultsTable);

%-------------------------------------------------------------------------------
% Compute a null distribution for each category size:
uniqueSizes = unique(resultsTable.size);
numSizes = length(uniqueSizes);
numGenes = height(geneInfo);

nullDistributionRaw = zeros(numSizes,numSamples);
nullDistributionAbs = zeros(numSizes,numSamples);
nullDistributionVE1 = zeros(numSizes,numSamples);
switch whatSurrogate
case 'independentSpatialShuffle'
    % Spatial shuffle case:
    fprintf(1,'(INDEPENDENT) SPATIAL SHUFFLE!!\n');
    params.g.humanOrMouse = sprintf('surrogate-%s',whatSpecies);
    parfor j = 1:numSamples
        fprintf(1,'Sample %u/%u\n',j,numSamples);
        % Get indpendently-shuffled data:
        [geneData,geneInfo,structInfo] = LoadMeG(params.g);
        for i = 1:numSizes
            shuffledData = geneData(:,1:uniqueSizes(i));
            [nullDistributionRaw(i,j),nullDistributionAbs(i,j),nullDistributionVE1(i,j)] = IntraCorrelationScore(shuffledData);
        end
    end
case 'geneShuffle'
    % Shuffle genes randomly:
    fprintf(1,'GENE SHUFFLE!!\n');
    parfor j = 1:numSamples
        fprintf(1,'Sample %u/%u\n',j,numSamples);
        % Get indpendently shuffled data:
        geneDataShuffle = geneData(:,randperm(numGenes));
        for i = 1:numSizes
            shuffledData = geneDataShuffle(:,1:uniqueSizes(i));
            [nullDistributionRaw(i,j),nullDistributionAbs(i,j),nullDistributionVE1(i,j)] = IntraCorrelationScore(shuffledData);
        end
    end
otherwise
    error('Unknown null model');
end

%-------------------------------------------------------------------------------
% Estimate p-values from a given test statistic:
switch pValsFromWhat
case 'raw'
    theField = sprintf('%s_raw',whatSpecies);
    theNullDistribution = nullDistributionAbs;
case 'abs'
    theField = sprintf('%s_abs',whatSpecies);
    theNullDistribution = nullDistributionAbs;
case 'VE1'
    theField = sprintf('%s_VE1',whatSpecies);
    theNullDistribution = nullDistributionVE1;
end

%-------------------------------------------------------------------------------
% Plot:
f = figure('color','w'); hold('on')
h1 = histogram(resultsTable.(theField),'normalization','pdf');
h2 = histogram(theNullDistribution(1,:),'normalization','pdf');
h3 = histogram(theNullDistribution(end,:),'normalization','pdf');
legend([h1,h2,h3],{'realData',sprintf('nullSize%u',uniqueSizes(1)),sprintf('nullSize%u',uniqueSizes(end))});

%-------------------------------------------------------------------------------
% Compute p-values (bigger scores are better):
pVal = nan(numGOCategories,1);
pValZ = nan(numGOCategories,1);
parfor i = 1:numGOCategories
    categoryScore = resultsTable.(theField)(i);
    if ~isnan(categoryScore)
        nullForSize = theNullDistribution(resultsTable.size(i)==uniqueSizes,:);
        pVal(i) = mean(categoryScore < nullForSize);
        % Gaussian approximation:
        pValZ(i) = 1 - normcdf(categoryScore,mean(nullForSize),std(nullForSize));
    end
end

%-------------------------------------------------------------------------------
% Update the GO table:
resultsTable.pVal = pVal;
resultsTable.pValZ = pValZ;
% FDR-corrected values:
resultsTable.pValCorr = mafdr(pVal,'BHFDR','true');
resultsTable.pValZCorr = mafdr(pValZ,'BHFDR','true');

%-------------------------------------------------------------------------------
% Sort:
resultsTable = sortrows(resultsTable,'pValZ','ascend');
% resultsTable = sortrows(resultsTable,{'pVal','pValCorr','mouse_abs'},{'ascend','ascend','descend'});

end
