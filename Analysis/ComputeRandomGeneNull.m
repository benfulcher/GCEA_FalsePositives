function ComputeRandomGeneNull(myGOID)

myGOID = 61001;

%-------------------------------------------------------------------------------
% Bits of parameters:
params = GiveMeDefaultParams('mouse');
whatSpecies = params.g.humanOrMouse;
numNullMaps = 1000;
numNullSamplesPerMap = 100;
whatNullType = 'randomMap';
whatCorr = 'Spearman';
aggregateHow = 'mean';
saveOut = true;
beVerbose = true;

%-------------------------------------------------------------------------------
% Get real data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
numGenes = height(geneInfo);
numAreas = height(structInfo);

%-------------------------------------------------------------------------------
% Get a generic GO Table:
GOTable = GiveMeGOData(params,geneInfo.entrez_id);
numAreas = height(structInfo);
numGOCategories = height(GOTable);
numGenesReal = height(geneInfo);

%-------------------------------------------------------------------------------
% Match genes for this category:
theCategoryIndex = find(GOTable.GOID==myGOID)
theGenesEntrez = GOTable.annotations{theCategoryIndex};
matchMe = find(ismember(geneInfo.entrez_id,theGenesEntrez));
geneDataCategory = geneData(:,matchMe);
numGenesCategory = length(matchMe);
fprintf(1,'%u/%u genes from this GO category have matching records in the expression data\n',...
                    length(matchMe),length(theGenesEntrez));

%-------------------------------------------------------------------------------
% Under random phenotypes:
nullMaps = rand(numAreas,numNullMaps);

nullScores = cell(numNullMaps,1);
realScores = zeros(numNullMaps,1);
for i = 1:numNullMaps
    fprintf(1,'Null map %u/%u\n',i,numNullMaps);
    myNullMap = nullMaps(:,i);
    % Compute the random-gene null distribution for this map:
    nullScores{i} = zeros(numNullSamplesPerMap,1);
    for j = 1:numNullSamplesPerMap+1
        scoresForThisMap = zeros(numGenesCategory,1);
        for k = 1:numGenesCategory
            if j==1
                expressionVector = geneDataCategory(:,k); % Real data
            else
                % random gene
                theRandomGeneIndex = randi(size(geneDataCategory,2));
                expressionVector = geneDataCategory(:,theRandomGeneIndex);
            end
            scoresForThisMap(j) = corr(myNullMap,expressionVector,'type',whatCorr,'rows','pairwise');
        end
        % Aggregate across genes in the category as a mean:
        if j==1
            realScores(i) = nanmean(scoresForThisMap);
        else
            nullScores{i}(j-1) = nanmean(scoresForThisMap);
        end
    end
end



end
