function GOTable = NodeShuffleEnrichment(whatEnrichment,whatShuffle,numNulls,structFilter,whatSpecies,params)
% Idea is to shuffle node properties across nodes, generating a null
% distribution for each category separately

if nargin < 1
    whatEnrichment = 'degree';
end
if nargin < 2
    whatShuffle = 'twoIsocortex'; % 'anatomy', 'anatomyTwo', 'twoIsocortex', 'all'
end
if nargin < 3
    numNulls = 200;
end
if nargin < 4 || isempty(structFilter)
    structFilter = 'all'; % 'isocortex','all'
end
if nargin < 5 || isempty(whatSpecies)
    whatSpecies = 'mouse';
    fprintf(1,'Mouse by default\n');
end
if nargin < 6
    params = GiveMeDefaultParams(whatSpecies);
end

%-------------------------------------------------------------------------------
% Ensure no gene expression normalization:
% Warn if normalization applied:
if ~strcmp(params.g.normalizationGene,'none')
    warning('Gene expression normalization across genes?? (Consider -> ''none'') :-O')
end
if ~strcmp(params.g.normalizationRegion,'none')
    warning('Gene expression normalization across regions?? (Consider ->''none'') :-O')
end

%-------------------------------------------------------------------------------
% Load data using default settings:
%-------------------------------------------------------------------------------
[A_bin,regionAcronyms,adjPVals] = GiveMeAdj(params.c.connectomeSource,params.c.pThreshold,true,...
                                    params.c.whatWeightMeasure,params.c.whatHemispheres,params.c.structFilter);
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
[A_bin,geneData,structInfo,keepStruct] = filterStructures(structFilter,structInfo,A_bin,geneData);
GOTable = GetFilteredGOData(params.e.whatSource,params.e.processFilter,...
                                params.e.sizeFilter,geneInfo.entrez_id);
numGenes = size(geneData,1);
numGOCategories = height(GOTable);

%-------------------------------------------------------------------------------
% Define gene scoring function
%-------------------------------------------------------------------------------
switch whatEnrichment
case 'degree'
    % Score based on correlations to degree
    k = sum(A_bin,1)' + sum(A_bin,2);
    score_fn = @(x) corr(k,x,'type','Spearman','rows','pairwise');
    fprintf(1,'Scoring genes as SPEARMAN correlation coefficients with degree\n');
    fprintf(1,'(Degree computed from %u connections across %u areas)\n',...
                                    sum(A_bin(:)),length(A_bin));
case 'isocortex'
    isCTX = ismember(structInfo.divisionLabel,'Isocortex');
    score_fn = @(x) -log10(ranksum(x(isCTX),x(~isCTX)));
end

%-------------------------------------------------------------------------------
% Define shuffle function
%-------------------------------------------------------------------------------
switch whatShuffle
case 'twoIsocortex'
    shuffle_fn = @()AnatomyShuffle(structInfo.divisionLabel,'twoIsocortex');
case 'anatomyTwo'
    shuffle_fn = @()AnatomyShuffle(structInfo.divisionLabel,'twoBroad');
case 'anatomyFive'
    shuffle_fn = @()AnatomyShuffle(structInfo.divisionLabel,'fiveByEye');
case 'all'
    shuffle_fn = @()randperm(numGenes);
otherwise
    error('Unknown shuffle setting ''%s''',whatShuffle);
end
fprintf(1,'Shuffling over %u brain areas using ''%s''\n',...
                    length(shuffle_fn()),whatShuffle);

%-------------------------------------------------------------------------------
% Assign scores to categories of genes
%-------------------------------------------------------------------------------
categoryScores = nan(numGOCategories,numNulls+1);
parfor n = 1:numNulls+1
    fprintf(1,'Null %u/%u\n',n,numNulls+1);
    if n == 1
        permVector = 1:size(geneData,1);
    else
        permVector = shuffle_fn();
    end

    gScore = zeros(height(geneInfo),1);
    for i = 1:height(geneInfo)
        gScore(i) = score_fn(geneData(permVector,i));
    end

    % Record mean scores for each category:
    for j = 1:numGOCategories
        matchMe = ismember(geneInfo.entrez_id,GOTable.annotations{j});
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
GOTable = EstimatePVals(categoryScores,whatTail,GOTable);
ix_GO = ListCategories(geneInfo,GOTable);
GOTable = GOTable(ix_GO,:);
categoryScores = categoryScores(ix_GO);

numSig = sum(GOTable.pValZCorr < params.e.enrichmentSigThresh);
fprintf(1,'%u significant categories at p_corr < %.2f\n',numSig,params.e.enrichmentSigThresh);
display(GOTable(1:numSig,:));

NullSummaryPlots(GOTable,categoryScores);
SpecificNullPlots(GOTable,categoryScores);

end
