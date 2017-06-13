% Extract edge-based measures to correlate with gene coexpression
%-------------------------------------------------------------------------------
% Ben Fulcher 23-2-2016

pThreshold = 0.05;

%-------------------------------------------------------------------------------
% Load data:
[C,G,CDM] = LoadAllData('All');

%-------------------------------------------------------------------------------
% Get gene coexpression:
Corrs.uncorrected = GiveMeRegionCorrs(G,'norm','energy',C,'none');
Corrs.corrected = GiveMeRegionCorrs(G,'norm','energy',C,'distance','expFitAll');

%-------------------------------------------------------------------------------
% Get the adjacency matrix:
A = logical(GiveMeAdj('Oh',pThreshold,true,'','right'));
PlotColorMatrix(A,C.RegionStruct,[],[],'',0);
colormap(gray)

%-------------------------------------------------------------------------------
% Remove nodes with 0 out-degree:
k_out = sum(A,2); % out degree
deadEndNodes = (k_out==0);
A = A(~deadEndNodes,~deadEndNodes);
Corrs.corrected = Corrs.corrected(~deadEndNodes,~deadEndNodes);
Corrs.uncorrected = Corrs.uncorrected(~deadEndNodes,~deadEndNodes);
warning('Removed %u nodes with 0 out-degree :/',sum(k_out==0))
k_out = sum(A,2); % out degree
k_in = sum(A,1)'; % in degree
k_tot = k_out + k_in; % total degree
numNodes = length(k_out);

%-------------------------------------------------------------------------------
% Compute edge measures on the adjacency matrix:
[edgeBet,G,signalTraffic,signalCount] = GiveMeAllEdgeMeasures(A);

%-------------------------------------------------------------------------------
% Compare different classes of links:
kHub = 45;
T = tabulate(k_tot);
bar(T(:,1),T(:,2));
isHub = (k_tot > 45);
richMask = false(size(A)); richMask(isHub,isHub) = A(isHub,isHub);
feedInMask = false(size(A)); feedInMask(~isHub,isHub) = A(~isHub,isHub);
feedOutMask = false(size(A)); feedOutMask(isHub,~isHub) = A(isHub,~isHub);
feederMask = false(size(A)); feederMask(isHub,~isHub) = A(isHub,~isHub); feederMask(~isHub,isHub) = A(~isHub,isHub);
peripheralMask = false(size(A)); peripheralMask(~isHub,~isHub) = A(~isHub,~isHub);
% Plot the masks:
imagesc(richMask+feederMask*2+peripheralMask*3)
colormap([0,0,0;BF_getcmap('redyellowblue',3,0)])

%-------------------------------------------------------------------------------
% Differences in edge betweenness:
% q = log10(edgeBet);
% q = log10(G);
% q = (signalTraffic);
q = signalCount;

dataCell = {q(richMask),q(feederMask),q(peripheralMask)};
extraParams = struct();
extraParams.theColors = BF_getcmap('set1',3,1);
BF_JitteredParallelScatter(dataCell,[],[],[],extraParams)
ax = gca;
ax.XTick = 1:3;
ax.XTickLabel = {'rich','feeder','peripheral'};

dataCell = {q(richMask),q(feedInMask),q(feedOutMask),q(peripheralMask)};
extraParams = struct();
extraParams.theColors = BF_getcmap('set1',4,1);
BF_JitteredParallelScatter(dataCell,[],[],[],extraParams)
ax = gca;
ax.XTick = 1:4;
ax.XTickLabel = {'rich','feedin','feedout','peripheral'};

%-------------------------------------------------------------------------------
f = figure('color','w');
subplot(231); imagesc(log10(edgeBet)); title('log10 Edge betweenness'); axis square
subplot(232); imagesc(log10(G)); title('log10 communicability'); axis square
subplot(233); imagesc(log10(signalTraffic)); title('log10 signal traffic'); axis square;
subplot(234); imagesc(log10(signalCount)); title('log10 signal count'); axis square;
Corrs_plot = abs(Corrs.uncorrected);
Corrs_plot(A==0) = 0;
subplot(235); imagesc(Corrs_plot); axis square; title('Gene coexpression')
colormap([0,0,0;BF_getcmap('bluegreen',9,0)])
% plot(edgeBet(A),signalTraffic(A),'.k')

%-------------------------------------------------------------------------------
doWhat = 'corrected';
f = figure('color','w');
numThresholds = 9;
subplot(221); hold on; BF_PlotQuantiles(log10(edgeBet(A)),Corrs.(doWhat)(A),numThresholds,0,0); xlabel('log10 edge betweenness'); axis square
subplot(222); hold on; BF_PlotQuantiles(log10(G(A)),Corrs.(doWhat)(A),numThresholds,0,0); xlabel('log10 communicability'); axis square
subplot(223); hold on; BF_PlotQuantiles(log10(signalTraffic(A)),Corrs.(doWhat)(A),numThresholds,0,0); xlabel('log10 signal traffic'); axis square;
subplot(224); hold on; BF_PlotQuantiles(signalCount(A),Corrs.(doWhat)(A),numThresholds,0,0); xlabel('Signal count'); axis square;
ylabel(doWhat)
