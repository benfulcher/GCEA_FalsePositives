function resultsTable = IntraCorrelationByCategory(params,whatSurrogate,numSamples,pValsFromWhat,doSave)
% Annotate intra-category coexpression for different gene-expression datasets
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check input parameters:
if nargin < 1 || isempty(params)
    params = GiveMeDefaultParams('mouse');
end
if ischar(params)
    % Can specify species name:
    params = GiveMeDefaultParams(params);
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
if nargin < 5
    doSave = true;
end
if doSave
    doPlot = false;
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
    params.g.humanOrMouse = sprintf('surrogate-%s',params.humanOrMouse);
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
    theField = sprintf('intracorr_raw');
    theNullDistribution = nullDistributionRaw;
case 'abs'
    theField = sprintf('intracorr_abs');
    theNullDistribution = nullDistributionAbs;
case 'VE1'
    theField = sprintf('intracorr_VE1');
    theNullDistribution = nullDistributionVE1;
end

%-------------------------------------------------------------------------------
% Plot:
if doPlot
    f = figure('color','w'); hold('on')
    h1 = histogram(resultsTable.(theField),'normalization','pdf');
    h2 = histogram(theNullDistribution(1,:),'normalization','pdf');
    h3 = histogram(theNullDistribution(end,:),'normalization','pdf');
    legend([h1,h2,h3],{'realData',sprintf('nullSize%u',uniqueSizes(1)),sprintf('nullSize%u',uniqueSizes(end))});
end

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

% FDR-correct:
resultsTable.pValCorr = mafdr(pVal,'BHFDR','true');
resultsTable.pValZCorr = mafdr(pValZ,'BHFDR','true');

%-------------------------------------------------------------------------------
% Sort:
resultsTable = sortrows(resultsTable,'pValZ','ascend');
% resultsTable = sortrows(resultsTable,{'pVal','pValCorr','mouse_abs'},{'ascend','ascend','descend'});

%-------------------------------------------------------------------------------
% Save out:
if doSave
    fileOut = fullfile('DataOutputs',sprintf('Intra_%s_%s_%s_%u.mat',...
                        params.humanOrMouse,whatSurrogate,pValsFromWhat,numSamples))
    save(fileOut,'resultsTable','params','numSamples');
end

end
