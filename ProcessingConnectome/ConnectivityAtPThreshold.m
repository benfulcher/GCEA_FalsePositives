% Adj is weights and PVals is p-values
% Simple function zeros large p-values
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-06-02
% ------------------------------------------------------------------------------

function Adj_th = ConnectivityAtPThreshold(Adj,PVals,p_th,replaceWithWhat)
    
    if nargin < 2
        error(['Require input of adjacency matrix and corresponding p-values for each' ...
                        ' connection weight)']);
    end
    if nargin < 3 || isempty(p_th)
        p_th = 0.05;
    end
    if nargin < 4 || isempty(replaceWithWhat)
        replaceWithWhat = 0; % Give high-p links a weight of 0
    end

    Adj_th = Adj; % make a copy
    
    Adj_th((PVals > p_th)) = replaceWithWhat; % NaN elements with high p-values

end