function edgeMeasure = GiveMeEdgeMeasure(whatMeasure,A_bin,A_wei,onlyOnEdges,A_p)
% Computes a single specified edge measure
%-------------------------------------------------------------------------------

if nargin < 4
    onlyOnEdges = true; % only measure edge metrics on edges that exist
end
if nargin < 5
    A_p = []; % p-values assigned to each edge (in the case of the Oh et al.
              % linear regression model)
end


switch whatMeasure
case 'distance'
    warning('Distances only work at the moment for whole-brain data')
    C = load('Mouse_Connectivity_Data.mat','Dist_Matrix');
    edgeData = C.Dist_Matrix{1,1};
    if onlyOnEdges
        edgeData(A_bin==0) = 0;
    else
        edgeData(tril(true(size(edgeData)))) = 0;
    end
case 'connections'
    edgeMeasure = double(A_bin);
    edgeMeasure(eye(size(edgeMeasure))) = NaN; % exclude diagonal
case 'connected'
    % pairwise binary measure of connectivity
    edgeMeasure = double(A_bin | A_bin');
    edgeMeasure(tril(true(size(edgeMeasure)))) = NaN;
% BINARY:
case 'ktot_ktot'
    kin = sum(A_bin,1)';
    kout = sum(A_bin,2);
    ktot = kin + kout;
    edgeMeasure = ktot*ktot';
case 'kin_kin'
    kin = sum(A_bin,1)';
    edgeMeasure = kin*kin';
case 'kout_kout'
    kout = sum(A_bin,2);
    edgeMeasure = kout*kout';
case 'kin_kout'
    kin = sum(A_bin,1)';
    kout = sum(A_bin,2);
    edgeMeasure = kin*kout';
case 'kout_kin'
    kin = sum(A_bin,1)';
    kout = sum(A_bin,2);
    edgeMeasure = kout*kin';
case 'bin_edgeBet'
    edgeMeasure = edge_betweenness_bin(A_bin);
case 'log10_bin_edgeBet'
    edgeMeasure = log10(edge_betweenness_bin(A_bin));
case 'bin_communicability'
    edgeMeasure = communicability(A_bin);
case 'signalTraffic'
    % Markov signal traffic measure:
    edgeMeasure = SimulateMarkov(A_bin);
case 'signalCount'
    % Simulate dissipative propagating message:
    edgeMeasure = SimulateWalker(A_bin,0.9,2000);
% WEIGHTED:
case 'weight'
    edgeMeasure = A_wei;
case 'log10_weight'
    edgeMeasure = log10(A_wei);
case 'pValue'
    if isempty(A_p)
        error('Did not give edge-based p-value info, you dummy');
    end
    edgeMeasure = A_p;
case 'wei_edgeBet'
    % weighted edge betweenness
    edgeMeasure = edge_betweenness_wei(A_wei);
case 'wei_communicability'
    % weighted communicability
    edgeMeasure = communicability(A_wei);
case 'stot_stot'
    sin = sum(A_wei,1)';
    sout = sum(A_wei,2);
    stot = sin + sout;
    edgeMeasure = stot*stot';
case 'sin_sin'
    sin = sum(A_wei,1)';
    edgeMeasure = sin*sin';
case 'sout_sout'
    sout = sum(A_wei,2);
    edgeMeasure = sout*sout';
case 'sin_sout'
    sin = sum(A_wei,1)';
    sout = sum(A_wei,2);
    edgeMeasure = sin*sout';
case 'sout_sin'
    sin = sum(A_wei,1)';
    sout = sum(A_wei,2);
    edgeMeasure = sout*kin';
otherwise
    error('Unknown edge measure: %s',whatMeasure);
end

%-------------------------------------------------------------------------------
if onlyOnEdges
    fprintf(1,'Values are only be assigned to edges that exist\n');
    edgeMeasure(A_bin==0) = 0;
end

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
% kin = sum(A_bin,1)';
% kout = sum(A_bin,2);
% ktot = kin + kout;
%
% edgeMeasures.ktot_ktot = ktot*ktot';
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
% if onlyOnEdges
%     fprintf(1,'Values can only be assigned to edges that exist\n');
%     measureNames = fieldnames(edgeMeasures);
%     numMeasures = length(measureNames);
%     for i = 1:numMeasures
%         edgeMeasures.(measureNames{i})(A_bin==0) = 0;
%     end
% end

end
