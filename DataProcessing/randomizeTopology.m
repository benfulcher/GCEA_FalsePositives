function A_rand = randomizeTopology(Adj,numIter)
% Ben Fulcher
%-------------------------------------------------------------------------------

if nargin < 2
    numIter = 100;
end

[A_rand,eff] = randmio_dir(Adj,numIter);


end
