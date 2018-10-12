
randomizeHow = 'shuffleStructAll';
humanOrMouse = 'mouse';
%-------------------------------------------------------------------------------

params = GiveMeDefaultParams(humanOrMouse);
% Compute edge measures (+nulls):
edgeData = GiveMeDistanceMatrix(params.humanOrMouse);
nullEdgeData = SpatialShuffleNull(edgeData,randomizeHow,params.gcc.numNulls);

% We're symmetric (Euclidean distances) so we should only look at upper triangles:
fprintf(1,'Taking upper triangles of symmetric distance-based data...\n');
isUpperDiag = triu(true(size(edgeData)),+1);
edgeData = edgeData(isUpperDiag);
nullEdgeData = cellfun(@(x)x(isUpperDiag),nullEdgeData);

%-------------------------------------------------------------------------------
% Do the enrichment at the level of individual genes:
[GOTable,geneScores] = EdgeEnrichment(edgeData,nullEdgeData,randomizeHow,params);
