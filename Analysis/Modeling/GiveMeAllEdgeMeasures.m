function edgeMeasures = GiveMeAllEdgeMeasures(A_bin,A_weight)
% Computes a range of edge-based traffic measures

%-------------------------------------------------------------------------------
% Make binary version:
% Abin = (A > 0);

%-------------------------------------------------------------------------------
% Edge betweenness:
edgeMeasures.edgeBet = edge_betweenness_bin(A_bin);

% Communicability:
edgeMeasures.G = communicability(A_bin);

% Markov signal traffic measure:
edgeMeasures.signalTraffic = SimulateMarkov(A_bin);

% Simulate dissipative propagating message:
edgeMeasures.signalCount = SimulateWalker(A_bin,0.9,2000);

% Degree combinations:
kin = sum(A_bin,1)';
kout = sum(A_bin,2);
ktot = kin + kout;

edgeMeasures.ktot_ktot = ktot*ktot';
edgeMeasures.kin_kin = kin*kin';
edgeMeasures.kout_kout = kout*kout';
edgeMeasures.kin_kout = kin*kout';
edgeMeasures.kout_kin = kout*kin';

%-------------------------------------------------------------------------------
% WEIGHTED
%-------------------------------------------------------------------------------

% The weight itself!
edgeMeasures.weight = A_weight;

% Strength combinations
sin = sum(A_weight,1)';
sout = sum(A_weight,2);
stot = sin + sout;

edgeMeasures.stot_stot = stot*stot';
edgeMeasures.sin_sin = sin*sin';
edgeMeasures.sout_sout = sout*sout';
edgeMeasures.sin_sout = sin*sout';
edgeMeasures.sout_sin = sout*kin';

end
