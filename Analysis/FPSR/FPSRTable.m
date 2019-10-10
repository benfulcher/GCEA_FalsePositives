
%-------------------------------------------------------------------------------
% DATA LOADING (DO ONCE)
%-------------------------------------------------------------------------------
numNullSamples = 10000;
% Load in the null data:
GOTableNullMouseRandom = SurrogateEnrichmentProcess('mouse',numNullSamples,'randomUniform','');
GOTableNullMouseAC = SurrogateEnrichmentProcess('mouse',numNullSamples,'spatialLag','');
GOTableNullHuman = SurrogateEnrichmentProcess('human',numNullSamples,'randomUniform','');
GOTableNullHumanAC = SurrogateEnrichmentProcess('human',numNullSamples,'spatialLag','');
% Now the reference data:
GOTableNullMouseRef = SurrogateEnrichmentProcess('mouse',numNullSamples,'randomUniform','independentSpatialShuffle');
GOTableNullHumanRef = SurrogateEnrichmentProcess('human',numNullSamples,'randomUniform','independentSpatialShuffle');
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Some simple stats:

propMouseRef0 = mean(GOTableNullMouseRef.sumUnderSig==0);
propHumanRef0 = mean(GOTableNullHumanRef.sumUnderSig==0);


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
% Add mouse reference:
[~,ia,ib] = intersect(GOTableCombined.GOID,GOTableNullMouseRef.GOID);
GOTableCombined = GOTableCombined(ia,:);
GOTableCombined.refMouse = GOTableNullMouseRef.sumUnderSig(ib);
% Add human reference:
[~,ia,ib] = intersect(GOTableCombined.GOID,GOTableNullHumanRef.GOID);
GOTableCombined = GOTableCombined(ia,:);
GOTableCombined.refHuman = GOTableNullHumanRef.sumUnderSig(ib);

%-------------------------------------------------------------------------------
% Sort
%-------------------------------------------------------------------------------
myScore = GOTableCombined.sumUnderSigMouse + GOTableCombined.sumUnderSigMouseAC ...
            + GOTableCombined.sumUnderSigHuman + GOTableCombined.sumUnderSigHumanAC;
[~,ix] = sort(myScore,'descend');
GOTableCombined = GOTableCombined(ix,:);
display(GOTableCombined(1:100,:))


%-------------------------------------------------------------------------------
% Exhibited an increase:
didIncreaseMouseSBPrand = mean(GOTableCombined.sumUnderSigMouse>GOTableCombined.refMouse);
didIncreaseHumanSBPrand = mean(GOTableCombined.sumUnderSigHuman>GOTableCombined.refHuman);
% foldChangeMouse = GOTableCombined.sumUnderSigMouseAC./GOTableCombined.refMouse;

meanIncreaseMouseSBPrandSBPAC = mean(GOTableCombined.sumUnderSigMouseAC-GOTableCombined.sumUnderSigMouse)/10000;
meanIncreaseHumanSBPrandSBPAC = mean(GOTableCombined.sumUnderSigHumanAC-GOTableCombined.sumUnderSigHuman)/10000;

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

theCategoryID = 31638;
rowID = GOTableCombined.GOID==theCategoryID;

display(GOTableCombined(rowID,:));
% GOTableNullMouseRandom(,:)
% GOTableNullMouseAC(GOTableNullMouseAC.GOID==theCategoryID,:)
% GOTableNullHuman(GOTableNullHuman.GOID==theCategoryID,:)
% GOTableNullHumanAC(GOTableNullHumanAC.GOID==theCategoryID,:)

TellMeLiteratureStory(theCategoryID)
