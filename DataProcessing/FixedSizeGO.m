function GOTable = FixedSizeGO(GOTable,sizeFix)
% Returns a GO annotation table in which all categories have the same number of annotations
% (annotations are a random sample up to the maximum, sizeFix)
%-------------------------------------------------------------------------------

% For reproducibility:
rng(0,'twister');

% Loop over categories:
numGOCategories = height(GOTable);
for i = 1:numGOCategories
    annotationsHere = GOTable.annotations{i};
    numAnnotations = length(annotationsHere);
    if numAnnotations < sizeFix
        GOTable.annotations{i} = [];
        % shouldn't really happen—-assuming the table has been min-size filtered
        % at sizeFix
    else
        GOTable.annotations{i} = annotationsHere(randperm(numAnnotations,sizeFix));
    end
end

end
