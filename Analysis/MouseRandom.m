% MouseRandom
%-------------------------------------------------------------------------------

%===============================================================================
% Within-category correlation
%===============================================================================
numSamples = 100000;
params = GiveMeDefaultParams('mouse');
mouseIntra = IntraCorrelationByCategory(params,'geneShuffle',numSamples);
fileOut = fullfile('DataOutputs','mouseIntra_geneShuffle_20k.mat')
save(fileOut,'mouseIntra','params','numSamples');

%
% results.mouseIntra.pValCorr = results.mouseIntra.pValZCorr;
%
% %-------------------------------------------------------------------------------
% % Are coexpression scores related to spatial correlation scores?
% justMouseResults = struct('mouseDistance',results.mouseDistance,'mouseIntra',results.mouseIntra);
% % e.g., for mouse:
% [rowVectorResults,GOTerms,allGOIDs,tableNames] = CombineTables(justMouseResults,'mouse',...
%     {'pValZ','pValZ'});
% % {'meanScore','mouse'});
% plot(rowVectorResults(1,:),rowVectorResults(2,:),'.k')
% xlabel(tableNames{1})
% ylabel(tableNames{2})
