function [ I_optim, cost_optim] = optimizeCoilCurrentFrobeniusCondition( ...
    mrxisetup, L_dict, I_init, s_min, s_max, s_interSteps, maxIter, verbose, debug, parallel )
% optimizes the coil currents in MRXI by minimizing the Frobenius condition number of
% the MRXI Leadfield Matrix

% INPUT
% L_dict ... nCoils*nSensors x nVox; Dictionary leadfield matrix containing the leadfield matrices
% of all coils in use with individual activation and concatenated below
% each other
% I_init ... nCoils x nActivations; contains the initial guess for the coil
% currents
% s_min (optional) ... minimum stepsize, Default: 1e-3
% s_max (optional) ... maximum stepsize, Default: 1
% s_interSteps (optional) ... number of logarithmically equidistant steps
% between s_min and s_max, Default: 25
% maxIter (optional) ... maximum number of iterations (per activation if
% optimizationType = 'consecutive'), Default: Inf 
% optimizationType (optional) ... has to be either 'consecutive' or 'simultaneous' (recommended);
% determines if the activations are optimized one after another or
% simultaneously, Default: 'simultaneous'
% verbose (optional) ... if false, no output is given from the function, Default: true
% parallel (optional) ... if true, the optimization is done in parallel for
% multiple activations, Default: false

% OUTPUT
% I_optim ... nCoils x nActivations; optimized coil current matrix
% cost_optim ... value of the cost function using optimized coil currents
% I_storage ... nCoils x nActivations x nIterations; all steps from I_init
% to I_optim
% cost_storage ... 1 x nIterations; all costs from I_init to I_optim

narginMax = 10;

if nargin < narginMax
    parallel = false;
elseif isempty(parallel)
    parallel = false;
end

if nargin < narginMax-1
    debug = false;
elseif isempty(debug)
    debug = false;
end

if nargin < narginMax-2
    verbose = true;
elseif isempty(verbose)
    verbose = true;
end

if nargin < narginMax-3
    maxIter = Inf;
elseif isempty(maxIter)
    maxIter = Inf;
end

if nargin < narginMax-4
    s_interSteps = 25;
elseif isempty(s_interSteps)
    s_interSteps = 25;
end

if nargin < narginMax-5
    s_max = 1;
elseif isempty(s_max)
    s_max = 1;
end

if nargin < narginMax-6
    s_min = 1e-3;
elseif isempty(s_min)
    s_min = 1e-3;
end

[nCoils, nActivations] = size(I_init);
L_dict = L_dict./normest(L_dict);
nSensors = round(size(L_dict,1)/nCoils);
nVox = size(L_dict,2);

% get vectorized lead-field matrix
Lv = vectorizeLeadfieldMatrix(L_dict, 'f', nSensors);

% get derivative of gram matrix wrt the ith current
Li = zeros(nSensors,nVox,nCoils);
idx = 1:nSensors;
for c = 1:nCoils
    Li(:,:,c) = L_dict(idx,:);
    idx = idx + nSensors;
end


% initialize current matrix and normalize each activation to unit length
I = I_init;
I = I./norm(I);
I_grad = zeros(nCoils,nActivations);
% R = mrxisetup.getCoilRadii;

iter = 0;

s_search = logspace(log10(s_min), log10(s_max), s_interSteps);
cost_temp = zeros(s_interSteps, 1);

%%%%%%%%%%%%% for testing
if debug
    if nCoils~=3 || nActivations~=1
        error('Debugging only allowed for nCoils=3 and nActivations=1');
    end
    pointsperdim = 100;
    % radial coordinates
    rho = linspace(0, pi, pointsperdim);
    phi = linspace(0, pi, pointsperdim);
    r = 1;
    [PHI, RHO] = meshgrid(rho,phi);
    X = r.*sin(RHO).*cos(PHI);
    Y = r.*sin(RHO).*sin(PHI);
    Z = r.*cos(RHO);

    testcost = zeros(pointsperdim);
    testcond = testcost;
    ii = 0;
    for ir = 1:pointsperdim
        for ip = 1:pointsperdim
            ii = ii + 1;
            I_curr = [X(ip,ir);Y(ip,ir);Z(ip,ir)];        
            testcost(ip,ir) = algocost( Lv, I_curr, nSensors );
            testcond(ip,ir) = cond( vectorizeLeadfieldMatrix(Lv*I_curr, 'b', nSensors) );
        end
    end


    % visualization
    figure(122); clf;
    X_plot = [X;-X];
    Y_plot = [Y;-Y];
    Z_plot = [Z;-Z];
    costpoints_plot = [testcond;testcond];
    surf(X_plot,Y_plot,Z_plot,costpoints_plot);
    xlabel I_1; ylabel I_2; zlabel I_3
    daspect([1 1 1]);
    shading interp;
    c = colorbar;
    c.Label.String = 'Spectral condition';
    hold on;
    scatter3(I(1),I(2),I(3),'filled','k');
    %
    figure(123); clf;
    X_plot = [X;-X];
    Y_plot = [Y;-Y];
    Z_plot = [Z;-Z];
    costpoints_plot = [testcost;testcost];
    surf(X_plot,Y_plot,Z_plot,costpoints_plot);
    xlabel I_1; ylabel I_2; zlabel I_3
    daspect([1 1 1]);
    shading interp;
    c = colorbar;
    c.Label.String = 'Frobenius condition';
    hold on;
    scatter3(I(1),I(2),I(3),'filled','k');
    %
    if isinf(maxIter) || maxIter > 1e3
        I_plot=zeros(3,200);
    else
        I_plot=zeros(3,maxIter);
    end
    I_plot(:,1) = I;
    plotIdx = 1;
    nArcpoints = 10;
    arc = zeros(3,size(I_plot,2)*nArcpoints);
end
%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
L_curr = vectorizeLeadfieldMatrix(Lv*I, 'b', nSensors);
cost = algocost( Lv, I, nSensors );
if verbose
    fprintf('Initial cost: %.5e\n', cost);
end
iter = iter+1;
while true
    if iter > maxIter
        if verbose
            fprintf('Hit maximum number of iterations %.0f. Terminating optimization.\n', maxIter);
        end
        break
    end
    
    % compute Frobenius norm of L
    froL = norm(L_curr,'fro');
    
    % compute Pseudo
    pL_curr = pinv(L_curr);
    fropL = norm(pL_curr,'fro');
    
    % compute derivative of Frobenius condition number
    if parallel
        parfor a = 1:nActivations
            actidx = (a-1)*nSensors+1:a*nSensors;
            I_grad(:,a) = updateActivationGradientParallel(actidx, L_curr, pL_curr, Li, fropL, froL);
        end
    else
        I_grad = updateActivationGradient(nSensors, L_curr, pL_curr, L_dict, fropL, froL);
    end
    
%     I_grad = I_grad./repmat(R,1,nActivations);
    residual = norm(I_grad);
    I_grad = pinv(I_grad'); % for Gauss-Newton steps (usually works better), comment out for gradient descent
    I_grad = I_grad./norm(I_grad);

   % determine stepsize with largest cost reduction within s_min and s_max
    % using s_interSteps with logarithmically equal steps
    for si = 1:s_interSteps
        cost_temp(si) = algocost( Lv,  I - s_search(si)*I_grad, nSensors );
        if si > 1
            if cost_temp(si) > cost_temp(si-1)
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
        if verbose
            fprintf('Smallest stepsize %.3e increases cost. Terminating optimization.\n', s_search(1));
        end
        break;
    % else update stepsize and costs
    else
        s = s_search(si-1);
        I_curr = I - s*I_grad;
        I = I_curr./norm(I_curr);
        L_curr = vectorizeLeadfieldMatrix(Lv*I, 'b', nSensors);
        iter = iter+1;
        cost_prev = cost;
        cost = cost_temp(si-1);
    end
    
    %%%%%%%%%%%%% for testing
    if debug
        plotIdx = plotIdx+1;
        I_plot(:,plotIdx) = I;
        scatter3(I(1),I(2),I(3),'r');
        % store gradient-descent path as arc
        arcTemp = [ linspace(I_plot(1,plotIdx-1), I_plot(1,plotIdx), nArcpoints);
                linspace(I_plot(2,plotIdx-1), I_plot(2,plotIdx), nArcpoints);
                linspace(I_plot(3,plotIdx-1), I_plot(3,plotIdx), nArcpoints)];
        arc(:, (plotIdx-2)*nArcpoints+1:(plotIdx-1)*nArcpoints) = normcols(arcTemp);
    end
    %%%%%%%%%%%%%

    % print cost if demanded
    if verbose
        fprintf('currentCost: %.5e, step: %.3e, iter.: %.0f, costRed.: %3.4f%%, runtime: %.1f, residualSteepness: %.3e\n', cost, s, iter-1, 100*(cost_prev-cost)/cost_prev,toc, residual);
    end
end

I_optim = I;
cost_optim = cost;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print final cost if demanded
if verbose
    fprintf('Final cost: %.5e, runtime: %.1f\n', cost_optim, toc);
end

%%%%%%%%%%%%% for testing
if debug
    scatter3(I_optim(1), I_optim(2), I_optim(3), 'filled', 'g');
    figure(120);
    scatter3(I_optim(1), I_optim(2), I_optim(3), 'filled', 'g');
    arc(:,(plotIdx-1)*nArcpoints+1:end) = [];
    plot3(arc(1,:),arc(2,:),arc(3,:),'r','linewidth',2)
    figure(121);
    plot3(arc(1,:),arc(2,:),arc(3,:),'r','linewidth',2)
end
%%%%%%%%%%%%%
end


function cost = algocost( Lv, I, nSensors, tol )
    L = vectorizeLeadfieldMatrix(Lv*I, 'b', nSensors);
    if nargin == 4
        L = L./normest(L);
        pL = pinv( L ,tol);
    else
        pL = pinv(L);
    end
    cost = norm(L,'fro') * norm(pL,'fro');
end


function Ia = updateActivationGradientParallel(actidx, L_curr, pL_curr, Li, fropL, froL)
    % see function "updateActivationGradient" below for comments
    [nSens, nVox, nCoils] = size(Li);
    nAct = size(L_curr,1)/nSens;
    term1 = sum(sum(bsxfun(@times, Li, L_curr(actidx,:)),1),2)*fropL/froL;
    R = reshape(Li, nSens, nVox*nCoils);
    S = pL_curr(:,actidx)*R;
    Sr = reshape(pL_curr'*S, nSens*nAct, nVox, nCoils);
    term2 = sum(sum(bsxfun(@times, Sr, pL_curr'),1),2)*froL/fropL;
    Ia = squeeze(term1-term2);
end

function I_grad = updateActivationGradient(nSensors, L_curr, pL_curr, L_dict, fropL, froL)

    % get constants
    nVox = size(L_curr,2);
    nCoils = size(L_dict,1)/nSensors;
    nActs = size(L_curr,1)/nSensors;

    % calculate traces of all dictionary leadfield matrices with all
    % current matrices (left term of the equation)
    Lv_dict = vectorizeLeadfieldMatrix(L_dict,'f',nSensors);
    term1 = Lv_dict'*vectorizeLeadfieldMatrix(L_curr,'f',nSensors)*fropL/froL;
    
    % calculate all square matrices from all dictionary leadfield matrices times current pseudoinverses 
    R = reshape(Lv_dict, nSensors, nVox*nCoils);
    pLv_curr = vectorizeLeadfieldMatrix(pL_curr','f',nSensors);
    
    % calculate S=pL'*pL_k*L_dict_i for all nActs*nCoils combinations
    S = reshape(pLv_curr, nSensors, nVox*nActs)'*R;
    
    % reshape S into nVox x nVox x nCoils*nActs matrix for further
    % calculation
    S = reshape(permute(reshape(S,[nVox nActs nVox nCoils]), [1 3 4 2]), [nVox nVox nCoils*nActs]);
    
    % compute right term of gradient
    term2 = sum(sum(bsxfun(@times, S, pL_curr*pL_curr'),1),2)*froL/fropL;
    
    % calculate gradient matrix
    I_grad = term1(:)-term2(:);
    I_grad = reshape(I_grad,nCoils,nActs);
end
