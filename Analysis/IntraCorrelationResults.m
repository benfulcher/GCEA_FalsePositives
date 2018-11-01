%-------------------------------------------------------------------------------
% Annotate intra-category correlations for different gene expression datasets
%-------------------------------------------------------------------------------
results = struct();

%-------------------------------------------------------------------------------
% Mouse:
params = GiveMeDefaultParams('mouse');
params.c.structFilter = 'all';
results.mouse_all = AnnotateIntraCorrelations(params,[],'intraMouseAll');

%-------------------------------------------------------------------------------
% Human:
params = GiveMeDefaultParams('human');
params.c.structFilter = 'cortex';
results.humanCTX = AnnotateIntraCorrelations(params,[],'intraHumanCTX');

%-------------------------------------------------------------------------------
% Combine:
[GOID,ia,ib] = intersect(results.mouse_all.GOID,results.humanCTX.GOID);
intraMouse = results.mouse_all.intraMouseAll(ia);
intraHuman = results.humanCTX.intraHumanCTX(ib);

%-------------------------------------------------------------------------------
% Histogram
f = figure('color','w');
hold on
hMouse = histogram(intraMouse);
hHuman = histogram(intraHuman);
xlabel('Mean intra-category correlation')
legend([hMouse,hHuman],'mouseAll','humanCTX')

%-------------------------------------------------------------------------------
% Are correlations between genes within categories related between mouse and human?
f = figure('color','w');
plot(intraMouse,intraHuman,'.k')
[r,p] = corr(intraMouse,intraHuman,'type','Spearman');

%-------------------------------------------------------------------------------
combinedScore = intraMouse + intraHuman;
[~,ix] = sort(combinedScore,'descend');
GOTogether = results.mouse_all(ia,:);
GOTogether.intraHumanCTX = intraHuman;
display(GOTogether(ix(1:50),:))

%-------------------------------------------------------------------------------
% Order categories by highest intra-class correlation:
GOTable.meanScore = categoryScores;
GOTable = sortrows(GOTable,'meanScore','descend');
[~,ix] = sort(categoryScores,'descend');
