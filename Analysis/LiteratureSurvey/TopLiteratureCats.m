function TopLiteratureCats(pValCorrThreshold)
% Which GO categories are most frequently reported in the literature?
% (across surveys of both mouse and human studies)
%-------------------------------------------------------------------------------
if nargin < 1
    pValCorrThreshold = 0.05;
end

LitTableMouse = MakeLiteratureTable('mouse',pValCorrThreshold);
LitTableHuman = MakeLiteratureTable('human',pValCorrThreshold);

allGOIDs = union(LitTableHuman.GOID,LitTableMouse.GOID);
numGOTot = length(allGOIDs);

numMouseHumanStudies = zeros(numGOTot,2);
theStudies = cell(numGOTot,2);
for i = 1:numGOTot
    theSpotMouse = (LitTableMouse.GOID==allGOIDs(i));
    if any(theSpotMouse)
        numMouseHumanStudies(i,1) = LitTableMouse.numStudies(theSpotMouse);
        theStudies{i,1} = BF_cat(LitTableMouse.studyList{theSpotMouse});
    end
    theSpotHuman = (LitTableHuman.GOID==allGOIDs(i));
    if any(theSpotHuman)
        numMouseHumanStudies(i,2) = LitTableHuman.numStudies(theSpotHuman);
        theStudies{i,2} = BF_cat(LitTableHuman.studyList{theSpotHuman});
    end
end

% Sort on total of significant annotations
[~,ix] = sort(sum(numMouseHumanStudies,2),'descend');
numMouseHumanStudies = numMouseHumanStudies(ix,:);
theStudies = theStudies(ix,:);
allGOIDs = allGOIDs(ix);

%-------------------------------------------------------------------------------
% Get some info about em:
params = GiveMeDefaultParams();
params.e.sizeFilter = [0,1e6];
GOTableGeneric = GiveMeGOData(params);

numToMatch = size(numMouseHumanStudies,1);
categoryNames = cell(numToMatch,1);
IDLabels = cell(numToMatch,1);
for i = 1:numToMatch
    matchInd = find(GOTableGeneric.GOID==allGOIDs(i));
    if ~isempty(matchInd)
        categoryNames{i} = GOTableGeneric.GOName{matchInd};
        IDLabels{i} = GOTableGeneric.GOIDlabel{matchInd};
    else
        categoryNames{i} = 'Unknown';
        IDLabels{i} = '';
    end
end

numMouseStudies = numMouseHumanStudies(:,1);
numHumanStudies = numMouseHumanStudies(:,2);
TogetherTable = table(categoryNames,IDLabels,numMouseStudies,numHumanStudies);

display(TogetherTable(1:50,:))

%-------------------------------------------------------------------------------
% Make a table to write out:
IDLabel = IDLabels;
CategoryName = categoryNames;
ID = allGOIDs;
NumberOfHumanStudies = numHumanStudies;
NumberOfMouseStudies = numMouseStudies;
List_HumanStudyLabels = theStudies(:,2);
List_MouseStudyLabels = theStudies(:,1);
T = table(CategoryName,IDLabel,ID,NumberOfHumanStudies,NumberOfMouseStudies,...
            List_HumanStudyLabels,List_MouseStudyLabels);
fileOut = fullfile('SupplementaryTables',sprintf('LiteratureAnnotations_p0%02u.csv',pValCorrThreshold*100));
writetable(T,fileOut,'Delimiter',',','QuoteStrings',true);
fprintf(1,'Saved literature significance results to %s\n',fileOut);

end
