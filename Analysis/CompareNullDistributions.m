% function CompareNullDistributions()

% load('connected-permutedGeneDep-biological_process-Gnone_Rnone-250nulls')
% load('wei_communicability-uniformTopology-biological_process-Gnone_Rnone-250nulls.mat');

% matFiles = {'connected-uniformTopology-biological_process-Gnone_Rnone-100nulls.mat',...
%         'connected-topology-biological_process-Gzscore_Rnone-100nulls.mat',...
%         'connected-permutedGeneDep-biological_process-Gnone_Rnone-250nulls.mat'};

matFiles = {'wei_communicability-permutedGeneDep-biological_process-Gnone_Rnone-250nulls.mat',...
            'wei_communicability-uniformTopology-biological_process-Gnone_Rnone-250nulls.mat'};

f = figure('color','w'); hold on

for i = 1:length(matFiles)
    load(matFiles{i},'meanNull');
    histogram(meanNull,'normalization','probability')
end
% Add the real mean:
load(matFiles{2},'categoryScores');
histogram(categoryScores(:,1),'normalization','probability')
% load('connected-permutedGeneDep-biological_process-Gnone_Rnone-250nulls.mat','categoryScores');
% histogram(categoryScores(:,1),'normalization','probability')
% plot(ones(2,1)*mean(categoryScores(:,1)),[0,max(get(gca,'YLim'))],'-r')
legend({'permutedGeneDep-mean','uniformTopology-mean','unpermuted'});

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
