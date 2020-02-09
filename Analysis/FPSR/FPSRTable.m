
%-------------------------------------------------------------------------------
% DATA LOADING (DO ONCE)
%-------------------------------------------------------------------------------
params = GiveMeDefaultParams('mouse');
numNullSamples = params.nulls.numNullsFPSR; % number of null maps to test against

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
%-------------------------------------------------------------------------------
% Proportion of the reference (proper null) data that have no false significance.
propMouseRef0 = mean(GOTableNullMouseRef.sumUnderSig==0);
propHumanRef0 = mean(GOTableNullHumanRef.sumUnderSig==0);
fprintf(1,'%u/%u (%g) mouse GO categories were never significant in the reference case\n',...
            sum(GOTableNullMouseRef.sumUnderSig==0),height(GOTableNullMouseRef),propMouseRef0);
fprintf(1,'%u/%u (%g) human GO categories were never significant in the reference case\n',...
            sum(GOTableNullHumanRef.sumUnderSig==0),height(GOTableNullHumanRef),propHumanRef0);

% Max FPSR:
maxRefMouse = max(GOTableNullMouseRef.sumUnderSig);
maxRefHuman = max(GOTableNullHumanRef.sumUnderSig);
fprintf(1,'Max FPSR (reference) of any mouse GO category was %g%%\n',...
                    maxRefMouse/numNullSamples*100);
fprintf(1,'Max FPSR (reference) of any human GO category was %g%%\n',...
                    maxRefHuman/numNullSamples*100);

% Mean FPSR:
meanRefMouse = mean(GOTableNullMouseRef.sumUnderSig);
meanRandMouse = mean(GOTableNullMouseRandom.sumUnderSig);
meanRefHuman = mean(GOTableNullHumanRef.sumUnderSig);
fprintf(1,'Mean FPSR (reference) of mouse GO categories was %g%%\n',...
                    meanRefMouse/numNullSamples*100);
fprintf(1,'Mean FPSR (SBP-random) of mouse GO categories was %g%%\n',...
                    meanRandMouse/numNullSamples*100);
fprintf(1,'Increase in mean FPSR (SBP-random) of mouse GO categories was %u-fold\n',...
                    round(meanRandMouse/meanRefMouse));
fprintf(1,'Mean FPSR (reference) of human GO categories was %g%%\n',...
                    meanRefHuman/numNullSamples*100);

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
% Save out to csv for paper:
IDLabel = GOTableCombined.GOIDlabel;
CategoryName = GOTableCombined.GOName;
ID = GOTableCombined.GOID;
FPSR_Mouse_Reference = GOTableCombined.refMouse/numNullSamples;
FPSR_Mouse_SBPrandom = GOTableCombined.sumUnderSigMouse/numNullSamples;
FPSR_Mouse_SBPspatial = GOTableCombined.sumUnderSigMouseAC/numNullSamples;
FPSR_Human_Reference = GOTableCombined.refHuman/numNullSamples;
FPSR_Human_SBPrandom = GOTableCombined.sumUnderSigHuman/numNullSamples;
FPSR_Human_SBPspatial = GOTableCombined.sumUnderSigHumanAC/numNullSamples;

T = table(CategoryName,IDLabel,ID,FPSR_Mouse_Reference,FPSR_Mouse_SBPrandom,...
                    FPSR_Mouse_SBPspatial,FPSR_Human_Reference,FPSR_Human_SBPrandom,...
                    FPSR_Human_SBPspatial);
fileOut = fullfile('SupplementaryTables','FPSRTable.csv');
writetable(T,fileOut,'Delimiter',',','QuoteStrings',true);
fprintf(1,'Saved all FPSR results to %s\n',fileOut);

%===============================================================================
% Some basic statistics on how FPSR estimates change across the scenarios
%-------------------------------------------------------------------------------
% Max:
maxMouseRand = max(GOTableCombined.sumUnderSigMouse);
maxHumanRand = max(GOTableCombined.sumUnderSigHuman);
maxMouseSpat = max(GOTableCombined.sumUnderSigMouseAC);
maxHumanSpat = max(GOTableCombined.sumUnderSigHumanAC);

% Exhibited an increase:
didIncreaseMouseSBPrand = mean(GOTableCombined.sumUnderSigMouse > GOTableCombined.refMouse);
didIncreaseHumanSBPrand = mean(GOTableCombined.sumUnderSigHuman > GOTableCombined.refHuman);
% foldChangeMouse = GOTableCombined.sumUnderSigMouseAC./GOTableCombined.refMouse;

meanIncreaseMouseSBPrandSBPAC = mean(GOTableCombined.sumUnderSigMouseAC-GOTableCombined.sumUnderSigMouse)/numNullSamples;
meanIncreaseHumanSBPrandSBPAC = mean(GOTableCombined.sumUnderSigHumanAC-GOTableCombined.sumUnderSigHuman)/numNullSamples;

%===============================================================================
% % Look up a specific category:
%
% %-------------------------------------------------------------------------------
% % ~~~~Do this once~~~~:
% % Get generic mouse-annotation GO category sizes:
% params = GiveMeDefaultParams();
% params.e.sizeFilter = [0,1e6];
% GOTerms = GiveMeGOData(params);
%
% %-------------------------------------------------------------------------------
% theCategoryName = 'long-term synaptic potentiation';
%
% weAreHere = strcmp(GOTerms.GOName,theCategoryName);
% display(GOTerms(weAreHere,:));
% theCategoryID = GOTerms.GOID(weAreHere);
%
% theCategoryID = 31638;
% rowID = GOTableCombined.GOID==theCategoryID;
%
% display(GOTableCombined(rowID,:));
% % GOTableNullMouseRandom(,:)
% % GOTableNullMouseAC(GOTableNullMouseAC.GOID==theCategoryID,:)
% % GOTableNullHuman(GOTableNullHuman.GOID==theCategoryID,:)
% % GOTableNullHumanAC(GOTableNullHumanAC.GOID==theCategoryID,:)
%
% TellMeLiteratureStory(theCategoryID)
