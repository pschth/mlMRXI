function M = optimizePreconditionerWRTSensitivity( L, s_max, s_min, s_interSteps, verbose )
% optimizes the coil positions using the Frobenius norm of the Gram matrix
% of the Leadfield matrix as cost function to be minimized.
% 
% INPUT:
% SETUP - MRXI setup configuration object
% 
% InitialPositions - 2 x nCoils matrix; initial coil positions specified by
% angles (radians) in the first matrix row and the distances along the vector from
% the bottom to the top of the cylinder
% 
% s_init - scalar; initial gradient-descent stepsize (recommendation: 1e-1)
% 
% s_thresh - scalar; lower stepsize threshold (recommendation: 1e-3)
% 
% verbose (optional) - boolean; prints optimization step information;
% default: true

nVox = size(L,2);

% randomly initialize preconditioning matrix
M = randn(size(L,1));
M = M./normest(M);

% initialize loop
grad = zeros(size(M));
cost = froGramCost( M*L ); % cost of current step
if verbose
    fprintf('Initial cost: %.6f, cost of unchanged L: %.6f\n', cost, froGramCost(L));
end

tic
while true
    % precalculate some stuff required for gradient calculation
    MtM = M'*M;
    litMtMli = zeros(nVox,1);
    for i = 1:nVox
        litMtMli(i) = L(:,i)'*MtM*L(:,i);
    end

    % calculate gradient
    for i = 1:nVox-1
        LitMtM = L(:,i)'*MtM;
        for j = i+1:nVox
            LiLjt = L(:,i)*L(:,j)';
            grad = grad + ...
                (LiLjt + LiLjt')*(LitMtM*L(:,j));
        end
    end
    grad = M*grad;
    grad = grad./normest(grad);
    
    % determine stepsize with largest cost reduction within s_min and s_max
    % using s_interSteps with logarithmically equal 
    s_search = logspace(log10(s_min), log10(s_max), s_interSteps);
    cost_temp = zeros(s_interSteps, 1);
    for si = 1:s_interSteps
        cost_temp(si) = froGramCost( (M - s_search(si)*grad) * L );
        if si > 1
            if cost_temp(si) > cost_temp(si-1)
                s_search(si:end) = [];
                cost_temp(si:end) = [];
                break;
            end
        else
            if cost_temp(1) > cost
                break;
            end
        end
    end
    
    % if no reduction of cost was found within the search space, exit loop
    if si == 1
        break;
    % else update stepsize and costs
    else
        s = s_search(end);
        cost = cost_temp(end);
    end
    
    % update preconditioner
    M = M - s*grad;
    
    % print cost if demanded
    if verbose
        fprintf('Current cost: %.6f, stepsize: %.3e, runtime: %.1f\n', cost, s, toc);
    end
    
end
% print final cost if demanded
if verbose
    fprintf('Final cost: %.6f, runtime: %.1f\n', cost, toc);
end



