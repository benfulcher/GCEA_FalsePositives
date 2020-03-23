
%-------------------------------------------------------------------------------
% DATA LOADING (DO ONCE)
%-------------------------------------------------------------------------------

% ############
% MOUSE:
% ############
params = GiveMeDefaultParams('mouse');
% Reference:
params.g.whatSurrogate = 'randomMap';
params.nulls.customShuffle = 'independentSpatialShuffle';
GOTableNullMouseRef = SurrogateEnrichmentProcess(params);
% SBP-random:
params.g.whatSurrogate = 'randomMap';
params.nulls.customShuffle = 'none';
GOTableNullMouseRandom = SurrogateEnrichmentProcess(params);
% SBP-spatial:
params.g.whatSurrogate = 'spatialLag';
params.nulls.customShuffle = 'none';
GOTableNullMouseAC = SurrogateEnrichmentProcess(params);

% ############
% HUMAN:
% ############
params = GiveMeDefaultParams('human');
% Reference:
params.g.whatSurrogate = 'randomMap';
params.nulls.customShuffle = 'independentSpatialShuffle';
GOTableNullHumanRef = SurrogateEnrichmentProcess(params);
% SBP-random:
params.g.whatSurrogate = 'randomMap';
params.nulls.customShuffle = 'none';
GOTableNullHumanRandom = SurrogateEnrichmentProcess(params);
% SBP-spatial:
params.g.whatSurrogate = 'spatialLag';
params.nulls.customShuffle = 'none';
GOTableNullHumanAC = SurrogateEnrichmentProcess(params);
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

% Max CFPR across categories:
maxRefMouse = max(GOTableNullMouseRef.sumUnderSig);
maxRefHuman = max(GOTableNullHumanRef.sumUnderSig);
fprintf(1,'Max CFPR (reference) of any mouse GO category is %g%%\n',...
                    maxRefMouse/numNullSamples*100);
fprintf(1,'Max CFPR (reference) of any human GO category is %g%%\n',...
                    maxRefHuman/numNullSamples*100);

% Changes in mean CFPR across all categories:
meanRefMouse = mean(GOTableNullMouseRef.sumUnderSig);
meanRandMouse = mean(GOTableNullMouseRandom.sumUnderSig);
meanACMouse = mean(GOTableNullMouseAC.sumUnderSig);
meanRefHuman = mean(GOTableNullHumanRef.sumUnderSig);
meanRandHuman = mean(GOTableNullHumanRandom.sumUnderSig);
meanACHuman = mean(GOTableNullHumanAC.sumUnderSig);
fprintf(1,'Mean CFPR (reference) of MOUSE GO categories is %g%%\n',...
                    meanRefMouse/numNullSamples*100);
fprintf(1,'Mean CFPR (SBP-random) of MOUSE GO categories is %g%%\n',...
                    meanRandMouse/numNullSamples*100);
fprintf(1,'Mean CFPR (SBP-spatial) of MOUSE GO categories is %g%%\n',...
                    meanACMouse/numNullSamples*100);
fprintf(1,'Increase in mean CFPR (SBP-random) of MOUSE GO categories is %u-fold\n',...
                    round(meanRandMouse/meanRefMouse));
fprintf(1,'Mean CFPR (reference) of HUMAN GO categories is %g%%\n',...
                    meanRefHuman/numNullSamples*100);
fprintf(1,'Mean CFPR (SBP-random) of HUMAN GO categories is %g%%\n',...
                    meanRandHuman/numNullSamples*100);
fprintf(1,'Mean CFPR (SBP-spatial) of HUMAN GO categories is %g%%\n',...
                    meanACHuman/numNullSamples*100);
fprintf(1,'Increase in mean CFPR (SBP-random) of HUMAN GO categories is %u-fold\n',...
                    round(meanRandHuman/meanRefHuman));


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
[~,ia,ib] = intersect(GOTableCombined.GOID,GOTableNullHumanRandom.GOID);
GOTableCombined = GOTableCombined(ia,:);
GOTableCombined.sumUnderSigHuman = GOTableNullHumanRandom.sumUnderSig(ib);
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
% Estimate literature significance
%-------------------------------------------------------------------------------
% (inefficient: thousands of duplicate loading, but we can beef through)
numRows = height(GOTableCombined);
mouseLiterature = zeros(numRows,1);
humanLiterature = zeros(numRows,1);
for i = 1:numRows
    [mouseLiterature(i),humanLiterature(i)] = LiteratureMatches(GOTableCombined.GOID(i),false);
    fprintf(1,'%u/%u\n',i,numRows);
end

%-------------------------------------------------------------------------------
% Save out to csv for paper:
IDLabel = GOTableCombined.GOIDlabel;
CategoryName = GOTableCombined.GOName;
ID = GOTableCombined.GOID;
CFPR_Mouse_Reference = GOTableCombined.refMouse/numNullSamples;
CFPR_Mouse_SBPrandom = GOTableCombined.sumUnderSigMouse/numNullSamples;
CFPR_Mouse_SBPspatial = GOTableCombined.sumUnderSigMouseAC/numNullSamples;
CFPR_Human_Reference = GOTableCombined.refHuman/numNullSamples;
CFPR_Human_SBPrandom = GOTableCombined.sumUnderSigHuman/numNullSamples;
CFPR_Human_SBPspatial = GOTableCombined.sumUnderSigHumanAC/numNullSamples;

T = table(CategoryName,IDLabel,ID,CFPR_Mouse_Reference,CFPR_Mouse_SBPrandom,...
                    CFPR_Mouse_SBPspatial,mouseLiterature,...
                    CFPR_Human_Reference,CFPR_Human_SBPrandom,...
                    CFPR_Human_SBPspatial,humanLiterature);
fileOut = fullfile('SupplementaryTables','CFPR_Table.csv');
writetable(T,fileOut,'Delimiter',',','QuoteStrings',true);
fprintf(1,'Saved all CFPR results to %s\n',fileOut);

%===============================================================================
% Some basic statistics on how CFPR estimates change across the scenarios
%-------------------------------------------------------------------------------
% Max:
maxMouseRand = max(GOTableCombined.sumUnderSigMouse);
maxHumanRand = max(GOTableCombined.sumUnderSigHuman);
maxMouseSpat = max(GOTableCombined.sumUnderSigMouseAC);
maxHumanSpat = max(GOTableCombined.sumUnderSigHumanAC);

% Exhibited an increase:
didIncreaseMouseSBPrand = mean(GOTableCombined.sumUnderSigMouse > GOTableCombined.refMouse);
didIncreaseHumanSBPrand = mean(GOTableCombined.sumUnderSigHuman > GOTableCombined.refHuman);
fprintf(1,'MOUSE: Ref -> SBP-rand, %.2f%% categories increased CFPR\n',didIncreaseMouseSBPrand*100);
fprintf(1,'HUMAN: Ref -> SBP-rand, %.2f%% categories increased CFPR\n',didIncreaseHumanSBPrand*100);
% foldChangeMouse = GOTableCombined.sumUnderSigMouseAC./GOTableCombined.refMouse;

meanIncreaseMouseSBPrandSBPAC = mean(GOTableCombined.sumUnderSigMouseAC-GOTableCombined.sumUnderSigMouse)/numNullSamples;
meanIncreaseHumanSBPrandSBPAC = mean(GOTableCombined.sumUnderSigHumanAC-GOTableCombined.sumUnderSigHuman)/numNullSamples;
fprintf(1,'MOUSE: Mean increase SBP-rand -> SBP-spatial = %.2f%%\n',meanIncreaseMouseSBPrandSBPAC*100);
fprintf(1,'HUMAN: Mean increase SBP-rand -> SBP-spatial = %.2f%%\n',meanIncreaseHumanSBPrandSBPAC*100);

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
