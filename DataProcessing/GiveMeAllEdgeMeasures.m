function edgeMeasures = GiveMeAllEdgeMeasures(A_bin,A_wei,onlyOnEdges,A_p)
% Computes a range of edge-based traffic measures

if nargin < 3
    onlyOnEdges = true; % only measure edge metrics on edges that exist
end
if nargin < 4
    A_p = []; % p-values assigned to each edge (in the case of the Oh et al.
              % linear regression model)
end
%-------------------------------------------------------------------------------

measureNames = {'kin_kin',...
                'kout_kout',...
                'kin_kout',...
                'kout_kin',...
                'ktot_ktot',... % product of total degrees
                'bin_edgeBet',... % binary edge betweenness
                'bin_communicability',... % communicability
                'signalTraffic',...
                'signalCount',... % simulate dissipative propagating message:
                'weight',...
                'wei_edgeBet',... % weighted edge betweenness
                'wei_communicability',... % weighted communicability
                'stot_stot',... % strength products
                'sin_sin',...
                'sout_sout',...
                'sin_sout',...
                'sout_sin',...
                };

% The p-value
if ~isempty(A_p)
    measureNames{end+1} = 'pValue';
end

%-------------------------------------------------------------------------------
numMeasures = length(measureNames);

% Compute all measures using GiveMeEdgeStat:
edgeMeasures = struct();
for i = 1:numMeasures
    edgeMeasures.(measureNames{i}) = GiveMeEdgeStat(measureNames{i},A_bin,A_wei,onlyOnEdges,A_p);
end

end
