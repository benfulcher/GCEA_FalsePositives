function nullEdgeData = SpatialShuffleNull(edgeData,randomizeHow,numNulls,structInfo)
% Shuffles edge data according to a given (spatial) randomization method
% Permute gene profiles assigned to regions (later)

%-------------------------------------------------------------------------------
% Check inputs:
if nargin < 2
    randomizeHow = 'shuffleStructAll';
end
if nargin < 3
    numNulls = 100;
end

%-------------------------------------------------------------------------------
% Define the permutation function (anatomically constrained?):
switch randomizeHow
case 'shuffleStructAll'
    fprintf(1,'Shuffling all structures uniformly...\n');
    N = size(edgeData,1);
    shuffle_fn = @()randperm(N);
case 'shuffleStructTwo'
    fprintf(1,'Shuffling within and outside of cerebral cortex separately...\n');
    shuffle_fn = @()AnatomyShuffle(structInfo.divisionLabel,'twoBroad',false);
case 'shuffleStructFive'
    fprintf(1,'Shuffling within five brain divisions...\n');
    shuffle_fn = @()AnatomyShuffle(structInfo.divisionLabel,'fiveByEye',false);
otherwise
    error('Unknown spatial randomization method: ''%s''',randomizeHow);
end

%-------------------------------------------------------------------------------
% Compute permuted versions of the edge data:
nullEdgeData = cell(numNulls,1);
for i = 1:numNulls
    % Compute the desired edge measure:
    rp = shuffle_fn();
    nullEdgeData{i} = edgeData(rp,rp);
end

end
