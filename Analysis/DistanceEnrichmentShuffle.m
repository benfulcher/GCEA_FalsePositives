
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

edgeData(isUpperDiag) = 0;
numNulls = length(nullEdgeData);
for i = 1:numNulls
    nullEdgeData{i}(isUpperDiag) = 0;
end

%-------------------------------------------------------------------------------
% Do the enrichment at the level of individual genes:
[GOTable,geneScores] = EdgeEnrichment(edgeData,nullEdgeData,randomizeHow,params);

%-------------------------------------------------------------------------------
% Save to mat file:
fileName = sprintf('%s-%s-%s-G%s_R%s-%unulls.mat','distance',randomizeHow,...
        params.e.processFilter,params.g.normalizationGene,params.g.normalizationRegion,numNulls);
save(fullfile('DataOutputs',fileName));
fprintf(1,'Saved shuffled distance results to %s\n',fileName);
