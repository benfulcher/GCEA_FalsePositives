function signalTraffic = SimulateMarkov(A,doPlot)
% Simulates a Markov model on the adjacency matrix, A
% Returns a measure of edge traffic, signalTraffic

if nargin < 2
    doPlot = 0;
end

k_out = sum(A,2); % out degree
k_in = sum(A,1)'; % in degree
k_tot = k_out + k_in; % total degree
numNodes = length(k_out);

% Transition probabilities, W
D = diag(k_out);
% W = (A'*inv(D))'; % this is equivalent to D\A
W = D\A; % this is equivalent

%-------------------------------------------------------------------------------
% Simulate a stationary distribution to check the eigenvalue solution
% f = figure('color','w');
% p0 = rand(1,numNodes);
% p0 = p0/sum(p0); % ensure sum of 1
% p0(1) = 1;
% bar(p0(ix)); xlim([0.5,numNodes+0.5])
% pause(0.2)
% p = p0;
% for i = 1:100
%     p = p*W;
%     bar(p(ix)); drawnow
%     title(i)
%     pause(0.01)
% end

% Can also compute as first eigenvalue:
[V,D] = eig(W');
lambda = V(:,1)/sum(V(:,1));
signalTraffic = diag(lambda)*A; % multiply each row by the stationary distribution

if doPlot
    %-------------------------------------------------------------------------------
    % Plot:
    [~,ix] = sort(k_tot,'descend');
    f = figure('color','w'); colormap(BF_getcmap('redyellowblue',6,0))
    subplot(311);bar([k_in(ix),k_out(ix)],'stacked');xlim([1,numNodes]);
    subplot(312);imagesc(A(ix,ix));
    subplot(313);imagesc(log10(W(ix,ix)));

    %-------------------------------------------------------------------------------
    % Plot:
    f = figure('color','w');
    subplot(311);plot(k_out,p,'ok'); xlabel('k-out'); ylabel('stationary p')
    subplot(312);plot(k_in,p,'ok'); xlabel('k-in'); ylabel('stationary p')
    subplot(313);plot(k_tot,p,'ok'); xlabel('k-tot'); ylabel('stationary p')
end

end
