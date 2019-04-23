function myMatrixNull = ShuffleMyMatrix(myMatrix,shuffleHow)

% Get dimensions:
[numRows,numCols] = size(myMatrix);

% Shuffle:
switch shuffleHow
case 'randomUniform'
    % Uniformly distributed numbers between 0 and 1
    myMatrixNull = rand([numRows,numCols]);

case 'independentRowShuffle'
    % Surrogate maps generated through (independent) random shuffling across brain areas
    % (should be consistent with random noise)
    myMatrixNull = myMatrix;
    for j = 1:numCols
        rp = randperm(numRows);
        myMatrixNull(:,j) = myMatrix(rp,j);
    end

case 'coordinatedRowShuffle'
    % Surrogate maps generated through random shuffling across brain areas
    rp = randperm(numRows);
    myMatrixNull = myMatrix(rp,:);

case 'coordinatedColumnShuffle'
    % Random shuffling of genes (randomizing association with metadata)
    rp = randperm(numCols);
    myMatrixNull = myMatrix(:,rp);

end
