function edgeMeasures = GiveMeNullEdgeMeasures(randomizeHow,whatEdgeMeasure,A_bin,A_wei,numNulls,onlyOnEdges,structInfo)
% Idea is to produce null edge values for a given null
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% (structInfo required for anatomically constrained nulls)
%-------------------------------------------------------------------------------

N = length(A_bin);

switch randomizeHow
case 'uniformTopology'
    % Randomize the topology to have a uniform probability across all possible edges
    fprintf(1,'Randomizing edges uniformly\n');
    numLinks = sum(A_bin(:));

    edgeMeasures = cell(numNulls+1,1);
    for i = 1:numNulls+1
        if i == 1
            A_rand = A_bin;
        else
            A_rand_vector = zeros(N*(N-1),1);
            rp = randperm(N*(N-1));
            A_rand_vector(rp(1:numLinks)) = 1;
            A_bin_i = squareform(A_rand_vector(1:end/2));
            A_bin_ii = squareform(A_rand_vector(end/2+1:end));
            upper = triu(true(size(A_bin)),+1);
            A_rand = zeros(size(A_bin));
            A_rand(upper) = A_bin_i(upper);
            lower = tril(true(size(A_bin)),-1);
            A_rand(lower) = A_bin_ii(lower);
        end

        weightVector = A_wei(A_wei > 0);
        A_wei_rand = A_rand;
        A_wei_rand(A_wei_rand > 0) = weightVector;

        % Compute the desired edge measure:
        edgeMeasures{i} = GiveMeEdgeMeasure(whatEdgeMeasure,A_rand,A_wei_rand,onlyOnEdges);
    end

case 'topology'
    numIter = 50; % (num randomizations per edge)
    fprintf(1,'Generating %u degree (+in-strength) topological nulls using %u randomizations per edge\n',...
                        numNulls,numIter);
    edgeMeasures = cell(numNulls+1,1);
    for i = 1:numNulls+1
        if i == 1
            A_rand = A_bin;
            A_wei_rand = A_wei;
        else
            A_rand = randmio_dir(A_bin,numIter);
            A_wei_rand = randmio_dir(A_wei,numIter);
        end
        % Compute the desired edge measure:
        edgeMeasures{i} = GiveMeEdgeMeasure(whatEdgeMeasure,A_rand,A_wei_rand,onlyOnEdges);
    end

case {'shuffleStructAll','shuffleStructTwo','shuffleStructFive'}
    % Permute gene profiles assigned to regions (later)
    % Edge measure stays constant using the correct topology:
    edgeMeasures0 = GiveMeEdgeMeasure(whatEdgeMeasure,A_bin,A_wei,onlyOnEdges);
    edgeMeasures = cell(numNulls+1,1);
    edgeMeasures{1} = edgeMeasures0;
    switch randomizeHow
    case 'shuffleStructAll'
        fprintf(1,'Shuffling all structures uniformly...\n');
        shuffle_fn = @()randperm(N);
    case 'shuffleStructTwo'
        fprintf(1,'Shuffling within and outside of cerebral cortex separately...\n');
        shuffle_fn = @()AnatomyShuffle(structInfo.divisionLabel,'twoBroad',false);
    case 'shuffleStructFive'
        fprintf(1,'Shuffling within five brain divisions...\n');
        shuffle_fn = @()AnatomyShuffle(structInfo.divisionLabel,'fiveByEye',false);
    end
    for i = 2:numNulls+1
        % Compute the desired edge measure:
        rp = shuffle_fn();
        edgeMeasures{i} = edgeMeasures0(rp,rp);
    end

case 'shuffleEdgeVals'
    % Shuffle values given to each edge
    % First compute the proper values:
    edgeMeasure0 = GiveMeEdgeMeasure(whatEdgeMeasure,A_bin,A_wei,onlyOnEdges);
    % Get the vector of values:
    edgeValsVector = edgeMeasure0(edgeMeasure0 > 0);
    % Now shuffle the values many times across the edges, preserving the binary topology
    edgeMeasures = cell(numNulls+1,1);
    edgeMeasures{1} = edgeMeasure0;
    for i = 2:numNulls+1
        edgeMeasures{i} = zeros(size(A_bin));
        edgeMeasures{i}(A_bin > 0) = edgeValsVector(randperm(length(edgeValsVector)));
    end
end
