%-------------------------------------------------------------------------------
% IntraCorrelationResults
%-------------------------------------------------------------------------------
% Annotate intra-category correlations for different gene expression datasets
%-------------------------------------------------------------------------------
numNullSamples_VE1 = 20000; % (Intra_*_VE1_20000.mat)
numNullSamples_surrogate = 10000; % (SurrogateGOTables_10000_*.mat)
whatShuffle = 'geneShuffle'; % 'geneShuffle', 'independentSpatialShuffle'


fileNameIn = @(whatSpecies) sprintf('Intra_%s_%s_VE1_%u.mat',whatSpecies,whatShuffle,numNullSamples_VE1);
%-------------------------------------------------------------------------------
results = struct();
resultsIntraMouse = load(fileNameIn('mouse'));
results.mouse = resultsIntraMouse.resultsTable;
resultsIntraHuman = load(fileNameIn('human'));
results.human = resultsIntraHuman.resultsTable;

[commonGOIDs,ia,ib] = intersect(results.mouse.GOID,results.human.GOID);

GOName = results.mouse.GOName(ia);
GOIDlabel = results.mouse.GOIDlabel(ia);
GOID = commonGOIDs;
sizeMouse = results.mouse.size(ia);
sizeHuman = results.human.size(ib);
intraVE1_mouse = results.mouse.intracorr_VE1(ia);
intraVE1_human = results.human.intracorr_VE1(ib);
newTable = table(GOName,GOIDlabel,GOID,sizeMouse,sizeHuman,intraVE1_mouse,intraVE1_human);
VE1_sum = intraVE1_mouse + intraVE1_human;
[~,ix] = sort(VE1_sum,'descend');
newTable = newTable(ix,:);

%===============================================================================
%===============================================================================
%===============================================================================
% OLDER STUFF:
%===============================================================================
%===============================================================================
%===============================================================================

%-------------------------------------------------------------------------------
% Get default parameters:
mouseParams = GiveMeDefaultParams('mouse');
humanParams = GiveMeDefaultParams('human');

%-------------------------------------------------------------------------------
% Mouse:
results.mouse_all = AnnotateIntraCorrelations(mouseParams,[],'mouse');

%-------------------------------------------------------------------------------
% Mouse nulls
params = mouseParams;
params.g.humanOrMouse = 'surrogate-mouse';

%       (i) random noise maps:
params.g.whatSurrogate = 'independentSpatialShuffle';
results.mouse_SpatialNull = AnnotateIntraCorrelations(params,[],'mouseSpatialShuffle');

%       (ii) gene metadata shuffled; genes assigned to categories at random:
params.g.whatSurrogate = 'geneShuffle';
results.mouse_geneNull = AnnotateIntraCorrelations(params,[],'mouseGeneShuffle');

% Check they all match:
if ~all(results.mouse_SpatialNull.GOID==results.mouse_all.GOID) || ...
        ~all(results.mouse_geneNull.GOID==results.mouse_all.GOID)
    error('Bad matching');
end

%-------------------------------------------------------------------------------
% Human:
results.humanCTX = AnnotateIntraCorrelations(humanParams,[],'human');

%-------------------------------------------------------------------------------
% Human nulls
params = humanParams;
params.g.humanOrMouse = 'surrogate-human';

%       (i) random noise maps:
params.g.whatSurrogate = 'independentSpatialShuffle';
results.human_SpatialNull = AnnotateIntraCorrelations(params,[],'humanSpatialShuffle');

%       (ii) gene metadata shuffled; genes assigned to categories at random:
params.g.whatSurrogate = 'geneShuffle';
results.human_geneNull = AnnotateIntraCorrelations(params,[],'humanGeneShuffle');

% Check they all match:
if ~all(results.human_SpatialNull.GOID==results.humanCTX.GOID) || ...
        ~all(results.human_geneNull.GOID==results.humanCTX.GOID)
    error('Bad matching');
end

%-------------------------------------------------------------------------------
% Combine:
[GOID,ia,ib] = intersect(results.mouse_all.GOID,results.humanCTX.GOID);
doAbs = true;
if doAbs
    intraMouse = results.mouse_all.mouse_abs(ia);
    intraMouseSpatialNull = results.mouse_SpatialNull.mouseSpatialShuffle_abs(ia);
    intraMouseGeneNull = results.mouse_geneNull.mouseGeneShuffle_abs(ia);
    intraHuman = results.humanCTX.human_abs(ib);
    intraHumanGeneNull = results.human_geneNull.humanGeneShuffle_abs(ib);
    intraHumanSpatialNull = results.human_SpatialNull.humanSpatialShuffle_abs(ib);
else
    intraMouse = results.mouse_all.mouse(ia);
    intraMouseSpatialNull = results.mouse_SpatialNull.mouseSpatialShuffle(ia);
    intraMouseGeneNull = results.mouse_geneNull.mouseGeneShuffle(ia);
    intraHuman = results.humanCTX.human(ib);
    intraHumanGeneNull = results.human_geneNull.humanGeneShuffle(ib);
    intraHumanSpatialNull = results.human_SpatialNull.humanSpatialShuffle(ib);
end

%-------------------------------------------------------------------------------
% Jittered scatter of distributions:
extraParams = struct();
extraParams.theColors = {[32,178,170]/255, [119,136,153]/255, [220,220,220]/255, [184,134,11]/255, [119,136,153]/255, [220,220,220]/255};
extraParams.customSpot = '';
extraParams.offsetRange = 0.7;
BF_JitteredParallelScatter({intraHuman,intraHumanGeneNull,intraHumanSpatialNull,...
                            intraMouse,intraMouseGeneNull,intraMouseSpatialNull},...
                            true,true,true,extraParams);
plot([0.5,6.5],[0,0],':k')
ax = gca; f = gcf();
ylabel('Mean intra-category correlation');
ax.XTick = 1:6;
ax.XTickLabel = {'human','human-geneShuffle','human-spatialShuffle','mouse','mouse-geneShuffle','mouse-spatialShuffle'};
f.Position = [1000,1099,305,239];
ax.XLim = [0.5,6.5];


% hHuman = histogram(intraHuman,'EdgeColor','k','FaceColor',/255);
% hMouse = histogram(intraMouse,'EdgeColor','k','FaceColor',[]/255);
% legend([hMouse,hHuman],'mouseAll','humanCTX')

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
