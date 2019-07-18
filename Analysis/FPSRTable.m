


%-------------------------------------------------------------------------------
% DATA LOADING (DO ONCE)
%-------------------------------------------------------------------------------
% Load in the null data:
GOTableNullMouseRandom = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','');
GOTableNullMouseAC = SurrogateEnrichmentProcess('mouse',10000,'spatialLag','');
GOTableNullHuman = SurrogateEnrichmentProcess('human',10000,'randomUniform','');
GOTableNullHumanAC = SurrogateEnrichmentProcess('human',10000,'spatialLag','');
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% COMBINE:
%-------------------------------------------------------------------------------
% Combine mice:
[~,ia,ib] = intersect(GOTableNullMouseRandom.GOID,GOTableNullMouseAC.GOID);
GOTableCombined = GOTableNullMouseRandom(ia,:);
GOTableCombined.sumUnderSigMouseAC = GOTableNullMouseAC.sumUnderSig(ib);
GOTableCombined.Properties.VariableNames{'sumUnderSig'} = 'sumUnderSigMouse';
deleteCol = strcmp(GOTableCombined.Properties.VariableNames,'pValCorr');
GOTableCombined(:,deleteCol) = [];

% Add human
[~,ia,ib] = intersect(GOTableCombined.GOID,GOTableNullHuman.GOID);
GOTableCombined = GOTableCombined(ia,:);
GOTableCombined.sumUnderSigHuman = GOTableNullHuman.sumUnderSig(ib);
% Add human AC:
[~,ia,ib] = intersect(GOTableCombined.GOID,GOTableNullHumanAC.GOID);
GOTableCombined = GOTableCombined(ia,:);
GOTableCombined.sumUnderSigHumanAC = GOTableNullHumanAC.sumUnderSig(ib);

%-------------------------------------------------------------------------------
% Sort
%-------------------------------------------------------------------------------
myScore = GOTableCombined.sumUnderSigMouse + GOTableCombined.sumUnderSigMouseAC ...
            + GOTableCombined.sumUnderSigHuman + GOTableCombined.sumUnderSigHumanAC;
[~,ix] = sort(myScore,'descend');
GOTableCombined = GOTableCombined(ix,:);
display(GOTableCombined(1:100,:))

%===============================================================================
% Look up a specific category:

%-------------------------------------------------------------------------------
% ~~~~Do this once~~~~:
% Get generic mouse-annotation GO category sizes:
params = GiveMeDefaultParams();
params.e.sizeFilter = [0,1e6];
GOTerms = GiveMeGOData(params);

%-------------------------------------------------------------------------------
theCategoryName = 'long-term synaptic potentiation';

weAreHere = strcmp(GOTerms.GOName,theCategoryName);
display(GOTerms(weAreHere,:));
theCategoryID = GOTerms.GOID(weAreHere);

theCategoryID = 7612;
rowID = GOTableCombined.GOID==theCategoryID;

display(GOTableCombined(rowID,:));
% GOTableNullMouseRandom(,:)
% GOTableNullMouseAC(GOTableNullMouseAC.GOID==theCategoryID,:)
% GOTableNullHuman(GOTableNullHuman.GOID==theCategoryID,:)
% GOTableNullHumanAC(GOTableNullHumanAC.GOID==theCategoryID,:)

TellMeLiteratureStory(theCategoryID)
