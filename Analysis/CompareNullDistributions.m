% function CompareNullDistributions()

% load('connected-permutedGeneDep-biological_process-Gnone_Rnone-250nulls')
load('connected-uniformTopology-biological_process-Gnone_Rnone-100nulls.mat')
load('wei_communicability-uniformTopology-biological_process-Gnone_Rnone-250nulls.mat');

A = ([gScores{2:end}]);
figure('color','w'); histogram(mean(A,2))
figure('color','w'); histogram(mean(A,1))
[~,ix] = sort(mean(A,2),'ascend');
[~,ix0] = sort(abs(mean(A,2)),'ascend');

figure('color','w');
for i = 1:25
    subplot(5,5,i);
    histogram(geneData(:,ix(i)));
    if i==1, title('largest scores'); end;
end
figure('color','w');
for i = 1:25
    subplot(5,5,i);
    histogram(geneData(:,ix0(i)));
    if i==1, title('scores nearest 0'); end;
end

%-------------------------------------------------------------------------------
% Synthetic nulls:
numNulls = 1000;
split = 63;
myT = zeros(numNulls,1);
g = geneData(:,7);
GCC = g*g';
V = GCC(triu(true(size(GCC)),+1)); % GCC_upper
numV = length(V);
for i = 1:numNulls
    rp = randperm(numV);
    % % U-test between connected and unconnected:

    % Sampling through difference partitions:
    % [p,~,stats] = ranksum(V(rp(1:split)),V(rp(split+1:end)));

    % Sampling with replacement from the same distribution:
    % [p,~,stats] = ranksum(V(randi(numV,split,1)),V(randi(numV,numV-split,1)));

    % Sampling independently from Gaussian distributions:
    [p,~,stats] = ranksum(randn(split,1),randn(numV-split,1));

    % Normalized Mann-Whitney U test (given the sample size may change across features)
    n1 = split;
    n2 = numV - split;

    normuStat = (stats.ranksum - n1*(n1+1)/2)/n1/n2; % normalized uStat
    myT(i) = normuStat;
    % [h,pVal,~,stats] = ttest2(V(rp(1:split)),V(rp(split+1:end)),'Vartype','unequal');
    % myT(i) = stats.tstat;
end
f = figure('color','w'); hold on
histogram(myT,'normalization','probability')
plot(mean(myT)*ones(2,1),[0,max(get(gca,'YLim'))],'r')
