function [ I_optim, cost_optim] = optimizeCoilCurrentSensitivity( mrxisetup, L_dict, I_init, s_min, s_max, s_interSteps, maxIter, verbose, debug, S_desired )
% optimizes the coil currents in MRXI by minimizing the squared sum of
% voxel sensitivities

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
% S_desired (optional) ... nVoxel vector; encodes the desired sensitvity
% distribution; Default: ones(nVoxel,1)

% OUTPUT
% I_optim ... nCoils x nActivations; optimized coil current matrix
% cost_optim ... value of the cost function using optimized coil currents
% I_storage ... nCoils x nActivations x nIterations; all steps from I_init
% to I_optim
% cost_storage ... 1 x nIterations; all costs from I_init to I_optim

narginMax = 10;

if nargin < narginMax
    S_desired = ones(mrxisetup.getNumberOfVoxels,1);
elseif isempty(S_desired)
    S_desired = ones(mrxisetup.getNumberOfVoxels,1);
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

% reshape dictionary leadfield matrix into nSensor x nVoxel x nCoil matrix
L_dict = permute(reshape(L_dict', nVox, nSensors, nCoils), [2 1 3]);

% initialize current matrix and normalize each activation to unit length
I = I_init;
I = I./norm(I);
R = mrxisetup.getCoilRadii;

iter = 0;

s_search = logspace(log10(s_min), log10(s_max), s_interSteps);
cost_temp = zeros(s_interSteps, 1);

%%%%%%%%%%%%% for testing
if debug
    if nCoils~=3 || nActivations~=1
        error('Debugging only allowed for nCoils=3 and nActivations=1');
    end
    pointsperdim = 25;
    % radial coordinates
    rho = linspace(0, pi, pointsperdim);
    phi = linspace(0, pi, pointsperdim);
    r = 1;
    [PHI, RHO] = meshgrid(rho,phi);
    X = r.*sin(RHO).*cos(PHI);
    Y = r.*sin(RHO).*sin(PHI);
    Z = r.*cos(RHO);

    testcost = zeros(pointsperdim);
    ii = 0;
    for ir = 1:pointsperdim
        for ip = 1:pointsperdim
            ii = ii + 1;
            I_curr = [X(ip,ir);Y(ip,ir);Z(ip,ir)];        
            testcost(ip,ir) = algocost( Lv, I_curr, nSensors, S_desired );
        end
    end


    % visualization
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
    c.Label.String = 'Sensitivity';
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
cost = algocost( Lv, I, nSensors, S_desired );
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
    
    % compute derivative of Frobenius condition number
    I_grad = updateGradient(nSensors, L_curr, L_dict, S_desired);
    
%     I_grad = I_grad./repmat(R,1,nActivations);
    I_grad = I_grad./norm(I_grad).*norm(I);

   % determine stepsize with largest cost reduction within s_min and s_max
    % using s_interSteps with logarithmically equal steps
    for si = 1:s_interSteps
        I_curr = I - s_search(si)*I_grad;
%         I_curr = I_curr./norm(I_curr);
        cost_temp(si) = algocost( Lv,  I_curr, nSensors, S_desired );
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
        I = I - s*I_grad;
%         I = I_curr./norm(I_curr);
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
%         I_prev = I + s*I_grad;
%         I_prev = I_prev./norm(I_prev);
%         I_quiver = cross(cross(I_prev,I_grad),I_prev);
%         I_quiver = I_quiver./norm(I_quiver)*0.2;
%         quiver3(I_prev(1), I_prev(2), I_prev(3), ...
%                 -I_quiver(1), -I_quiver(2), -I_quiver(3), 'r');
        % store gradient-descent path as arc
        arcTemp = [ linspace(I_plot(1,plotIdx-1), I_plot(1,plotIdx), nArcpoints);
                linspace(I_plot(2,plotIdx-1), I_plot(2,plotIdx), nArcpoints);
                linspace(I_plot(3,plotIdx-1), I_plot(3,plotIdx), nArcpoints)];
        arc(:, (plotIdx-2)*nArcpoints+1:(plotIdx-1)*nArcpoints) = normcols(arcTemp);
    end
    %%%%%%%%%%%%%

    % print cost if demanded
    if verbose
        fprintf('Current cost: %.5e, stepsize: %.3e, iteration: %.0f, cost reduction: %3.4f%%, runtime: %.1f\n', cost, s, iter-1, 100*(cost_prev-cost)/cost_prev,toc);
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
    arc(:,(plotIdx-1)*nArcpoints+1:end) = [];
    plot3(arc(1,:),arc(2,:),arc(3,:),'r','linewidth',2)
end
%%%%%%%%%%%%%
end


function cost = algocost( Lv, I, nSensors, S_desired )
    % calculates the frame potential of leadfield matrix
    L = vectorizeLeadfieldMatrix(Lv*I, 'b', nSensors);
    cost = sum((sum(abs(L),1)'-S_desired).^2);
end

function I_grad = updateGradient(nSensors, L_curr, L_dict, S_desired)

    % get constants
    nVox = size(L_curr,2);
    nCoils = size(L_dict,3);
    nActs = size(L_curr,1)/nSensors;
    
    % get sensitivities of current leadfield matrix
    S = sum(abs(L_curr),1)';
    
    % reshape current leadfield matrix to nSensor x nVoxel x nActs matrix
    L_curr = permute(reshape(L_curr', nVox, nSensors, nActs), [2 1 3]);
    
    % initialize gradient
    I_grad = zeros(nCoils, nActs);
    
    % calculate partial gradients
    for a = 1:nActs
        for c = 1:nCoils
            I_grad(c,a) = sum(sign(L_curr(:,:,a)).*L_dict(:,:,c), 1) * (S-S_desired);
        end
    end
end
