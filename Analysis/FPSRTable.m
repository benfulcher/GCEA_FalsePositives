


%-------------------------------------------------------------------------------
% DATA LOADING (DO ONCE)
%-------------------------------------------------------------------------------
% Get generic mouse-annotation GO category sizes:
params = GiveMeDefaultParams();
params.e.sizeFilter = [0,1e6];
GOTerms = GiveMeGOData(params);
% Load in the null data:
GOTableNullMouseRandom = SurrogateEnrichmentProcess('mouse',10000,'randomUniform','')
GOTableNullMouseAC = SurrogateEnrichmentProcess('mouse',10000,'spatialLag','')
GOTableNullHuman = SurrogateEnrichmentProcess('human',10000,'randomUniform','')
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

theCategoryName = 'chemical synaptic transmission, postsynaptic';

weAreHere = strcmp(GOTerms.GOName,theCategoryName);
display(GOTerms(weAreHere,:));
theCategoryID = GOTerms.GOID(weAreHere);

GOTableNullMouseRandom(GOTableNullMouseRandom.GOID==theCategoryID,:)
GOTableNullMouseAC(GOTableNullMouseAC.GOID==theCategoryID,:)
GOTableNullHuman(GOTableNullHuman.GOID==theCategoryID,:)

LiteratureMatches(theCategoryID)
