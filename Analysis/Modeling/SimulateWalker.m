function signalCount = SimulateWalker(A,pp,numSims,beVerbose)
% Given an adjacency matrix, A, simulates a diffusion process on it
%-------------------------------------------------------------------------------

if nargin < 2
    pp = 0.9;
end
if nargin < 3
    numSims = 100;
end
if nargin < 4
    beVerbose = 0;
end
% fprintf(1,'Signal propagates with probability p = %.2f\n',p);

%-------------------------------------------------------------------------------
numNodes = length(A);
signalCount = zeros(size(A));

% Set up figure:
if beVerbose
    f = figure('color','w');
    colormap([0,0,0;flipud(BF_getcmap('redyellowblue',8,0))])
end

% Run iterations:
for i = 1:numSims
    fprintf(1,'Iteration %u/%u\n',i,numSims);

    % Generate a signal at a random node:
    hasSignal = zeros(1,numNodes);
    hasSignal(randi(numNodes)) = 1;

    p = 0.2; % initially, this probability for neighbors of a target to receive the signal
    r = 0; % counts the number of iterations the simulation continues for

    while any(hasSignal)
        % Propagate to out-neighbors with probability p (increase for multiple signals)

        signalOut = zeros(size(A));
        signalOut(hasSignal > 0,:) = diag(hasSignal(hasSignal>0))*A(hasSignal>0,:); % send multiple signals out if contain multiple signals
        % signalOut(hasSignal,:) = A(hasSignal,:); % max to 1 signal input
        totalOutputs = sum(signalOut(:) > 0);

        if beVerbose
            subplot(121); imagesc([repmat(hasSignal',1,10),signalOut]); axis square; colorbar
            title(sprintf('Output targets. Simulation %u (iteration %u), p = %.3f',i,r,p))
        end
        % Message is transferred with probability p:
        signalOut(signalOut > 0) = (rand(sum(signalOut(:) > 0),1) < (1-(1-p).^(signalOut(signalOut > 0))));

        % Tell the user about it:
        if beVerbose
            fprintf(1,'Propagated from %u active nodes (mean messages: %.2f) to %u nodes (/%u) (p = %.3f)\n',...
                        sum(hasSignal>0),mean(hasSignal(hasSignal>0)),sum(sum(signalOut,1) > 0),totalOutputs,p);
            subplot(122); imagesc([repmat(hasSignal',1,10),signalOut]); axis square; colorbar
            title('Links through which signal propagated')
            drawnow;
        end

        % Add new signal to total:
        signalCount = signalCount + signalOut;

        % Update hasSignal binary variable:
        hasSignal = sum(signalOut,1); % count how many signals input to this node

        % Reduce propagation probability (exponential decay until the system dies down):
        r = r+1;
        p = p*pp;
    end
end

% Normalize to mean across simulations:
signalCount = signalCount/numSims;

end
