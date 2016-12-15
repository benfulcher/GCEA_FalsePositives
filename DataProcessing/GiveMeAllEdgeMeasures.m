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
% Make binary version:
% Abin = (A > 0);

%-------------------------------------------------------------------------------
% % Edge betweenness:
% edgeMeasures.bin_edgeBet = edge_betweenness_bin(A_bin);
%
% % Communicability (assigned to all pairs of nodes):
% edgeMeasures.bin_communicability = communicability(A_bin);

% Markov signal traffic measure:
% edgeMeasures.signalTraffic = SimulateMarkov(A_bin);

% Simulate dissipative propagating message:
% edgeMeasures.signalCount = SimulateWalker(A_bin,0.9,2000);

% Degree combinations (assigned to all pairs of nodes):
kin = sum(A_bin,1)';
kout = sum(A_bin,2);
ktot = kin + kout;

edgeMeasures.ktot_ktot = ktot*ktot';
% edgeMeasures.kin_kin = kin*kin';
% edgeMeasures.kout_kout = kout*kout';
% edgeMeasures.kin_kout = kin*kout';
% edgeMeasures.kout_kin = kout*kin';

%-------------------------------------------------------------------------------
% WEIGHTED
%-------------------------------------------------------------------------------

% % The weight itself!
% edgeMeasures.weight = A_wei;
%
% % The p-value
% if ~isempty(A_p)
%     edgeMeasures.pValue = A_wei;
% end
%
% % Weighted edge betweenness
% edgeMeasures.wei_betweenness = edge_betweenness_wei(A_wei);
%
% % Weighted communicability:
% edgeMeasures.wei_communicability = communicability(A_wei);
%
% % Strength combinations (assigned to all pairs of nodes)
% sin = sum(A_wei,1)';
% sout = sum(A_wei,2);
% stot = sin + sout;
%
% edgeMeasures.stot_stot = stot*stot';
% edgeMeasures.sin_sin = sin*sin';
% edgeMeasures.sout_sout = sout*sout';
% edgeMeasures.sin_sout = sin*sout';
% edgeMeasures.sout_sin = sout*kin';

%-------------------------------------------------------------------------------
if onlyOnEdges
    fprintf(1,'Values can only be assigned to edges that exist\n');
    measureNames = fieldnames(edgeMeasures);
    numMeasures = length(measureNames);
    for i = 1:numMeasures
        edgeMeasures.(measureNames{i})(A_bin==0) = 0;
    end
end

end
