
randomizeHow = 'shuffleStructAll';
numNulls = 100;

% Compute edge measures (+nulls):
edgeData = GiveMeDistanceMatrix('mouse');
nullEdgeData = SpatialShuffleNull(edgeData,randomizeHow,numNulls);

% Gene expression data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
numGenes = height(geneInfo);
fprintf(1,'%u x %u gene expression matrix\n',size(geneData,1),size(geneData,2));

%
[GOTable,gScore] = EdgeEnrichment(edgeData,nullEdgeData,whatNull,params);
