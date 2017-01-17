% function CompareNullDistributions()

% load('connected-permutedGeneDep-biological_process-Gnone_Rnone-250nulls')
% load('wei_communicability-uniformTopology-biological_process-Gnone_Rnone-250nulls.mat');

matFiles = {'connected-uniformTopology-biological_process-Gnone_Rnone-100nulls.mat',...
        'connected-topology-biological_process-Gzscore_Rnone-100nulls.mat',...
        'connected-permutedGeneDep-biological_process-Gnone_Rnone-250nulls.mat'};

% load('connected-uniformTopology-biological_process-Gnone_Rnone-500nulls.mat','meanNull')
% load('connected-topology-biological_process-Gzscore_Rnone-100nulls.mat','meanNull')
% load('connected-permutedGeneDep-biological_process-Gnone_Rnone-2000nulls.mat','meanNull');

f = figure('color','w'); hold on

for i = 1:length(matFiles)
    load(matFiles{i},'meanNull');
    histogram(meanNull,'normalization','probability')
end
% Add the real mean:
load('connected-topology-biological_process-Gzscore_Rnone-100nulls.mat','categoryScores');
histogram(categoryScores(:,1),'normalization','probability')
load('connected-permutedGeneDep-biological_process-Gnone_Rnone-250nulls.mat','categoryScores');
histogram(categoryScores(:,1),'normalization','probability')
% plot(ones(2,1)*mean(categoryScores(:,1)),[0,max(get(gca,'YLim'))],'-r')
legend({'uniformTopology-mean','topology-mean','permutedGeneDep-mean',...
            'unpermuted category scores (1)','unpermuted category scores (2)'});

%-------------------------------------------------------------------------------

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
whatNull = 'randomTopology';
split = 63;
myTestStat = zeros(numNulls,1);
% g = geneData(:,7);
% g = BF_NormalizeMatrix(geneData(:,7),'mixedSigmoid');
g = BF_NormalizeMatrix(geneData(:,7),'zscore');
% g = randn(213,1);
% g = rand(213,1);
GCC = g*g';
V = GCC(triu(true(size(GCC)),+1)); % GCC_upper
% histogram(V)
numV = length(V);
rp0 = randperm(numV);
for i = 1:numNulls
    switch whatNull
    case 'randomTopology'
        % Random partition:
        rp = randperm(numV);
        %---------------------------------------------------------------------------
        % % U-test between connected and unconnected:

        % Sampling through difference partitions:
        % [p,~,stats] = ranksum(V(rp(1:split)),V(rp(split+1:end)));

        % Sampling with replacement from the same distribution:
        % [p,~,stats] = ranksum(V(randi(numV,split,1)),V(randi(numV,numV-split,1)));

        % Sampling independently from Gaussian distributions:
        % [p,~,stats] = ranksum(randn(split,1),randn(numV-split,1));

        % Normalized Mann-Whitney U test (given the sample size may change across features)
        % n1 = split;
        % n2 = numV - split;

        % normuStat = (stats.ranksum - n1*(n1+1)/2)/n1/n2; % normalized uStat
        % myT(i) = normuStat;

        %---------------------------------------------------------------------------
        % T_TEST
        X = V(rp(1:split));
        Y = V(rp(split+1:end));
        % [h,pVal,~,stats] = ttest2(X,Y,'Vartype','unequal');
        % myT(i) = stats.tstat;
    case 'permutedGeneDep'
        gt = g(randperm(size(g)));
        GCCt = g*g';
        V = GCCt(triu(true(size(GCCt)),+1)); % GCC_upper
        X = V(rp0(1:split)); % same every time
        Y = V(rp0(split+1:end)); % same permutation every time
    end
    myTestStat(i) = (nanmean(X)-nanmean(Y))/sqrt(nanstd(X)^2+nanstd(Y)^2);
end
f = figure('color','w'); hold on
histogram(myTestStat,'normalization','probability')
plot(mean(myTestStat)*ones(2,1),[0,max(get(gca,'YLim'))],'r')
