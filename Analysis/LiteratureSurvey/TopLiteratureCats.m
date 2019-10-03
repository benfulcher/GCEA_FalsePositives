%-------------------------------------------------------------------------------
% Which GO categories are most frequently reported in the literature?
%-------------------------------------------------------------------------------

pValCorrThreshold = 0.05;
LitTableHuman = MakeLiteratureTable('human',pValCorrThreshold);
LitTableMouse = MakeLiteratureTable('mouse',pValCorrThreshold);

allGOIDs = union(LitTableHuman.GOID,LitTableMouse.GOID);
numGOTot = length(allGOIDs);

numMouseHumanStudies = zeros(numGOTot,2);
for i = 1:numGOTot
    if any(LitTableMouse.GOID==allGOIDs(i))
        numMouseHumanStudies(i,1) = LitTableMouse.numStudies(LitTableMouse.GOID==allGOIDs(i));
    end
    if any(LitTableHuman.GOID==allGOIDs(i))
        numMouseHumanStudies(i,2) = LitTableHuman.numStudies(LitTableHuman.GOID==allGOIDs(i));
    end
end

% Sort on total
[~,ix] = sort(sum(numMouseHumanStudies,2),'descend');
numMouseHumanStudies = numMouseHumanStudies(ix,:);
allGOIDs = allGOIDs(ix);

%-------------------------------------------------------------------------------
% Get some info about em:
params = GiveMeDefaultParams();
params.e.sizeFilter = [0,1e6];
GOTableGeneric = GiveMeGOData(params);

numToMatch = size(numMouseHumanStudies,1);
categoryNames = cell(numToMatch,1);
for i = 1:numToMatch
    matchInd = find(GOTableGeneric.GOID==allGOIDs(i));
    if ~isempty(matchInd)
        categoryNames{i} = GOTableGeneric.GOName{matchInd};
    end
end

numMouseStudies = numMouseHumanStudies(:,1);
numHumanStudies = numMouseHumanStudies(:,2);
TogetherTable = table(categoryNames,numMouseStudies,numHumanStudies);

display(TogetherTable(1:50,:))
