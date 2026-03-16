function [ mrxisetup_optim, cost_optim] = optimizeCoilRadiiFrobeniusConditionWithLowerLimit(...
    mrxisetup_init, s_min, s_max, s_interSteps, maxIter, normL_flag, lowerRadiusLimit, verbose, debug )
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
% normL_flag (optional) ... boolean, if true, the individual system
% matrices of each coil are normalized such that Lc = Lc./norm(Lc,'fro').
% Default: false
% optimizationType = 'consecutive'), Default: Inf 
% optimizationType (optional) ... has to be either 'consecutive' or 'simultaneous' (recommended);
% determines if the activations are optimized one after another or
% simultaneously, Default: 'simultaneous'
% verbose (optional) ... if false, no output is given from the function, Default: true

% OUTPUT
% I_optim ... nCoils x nActivations; optimized coil current matrix
% cost_optim ... value of the cost function using optimized coil currents
% I_storage ... nCoils x nActivations x nIterations; all steps from I_init
% to I_optim
% cost_storage ... 1 x nIterations; all costs from I_init to I_optim

mrxisetup = mrxisetup_init.copySetup;

if nargin < 9
    debug = false;
elseif isempty(debug)
    debug = false;
end

if nargin < 8
    verbose = true;
elseif isempty(verbose)
    verbose = true;
end

if nargin < 7
    lowerRadiusLimit = -Inf;
elseif isempty(lowerRadiusLimit)
    lowerRadiusLimit = -Inf;
end

if nargin < 6
    normL_flag = false;
elseif isempty(normL_flag)
    normL_flag = false;
end

if nargin < 5
    maxIter = Inf;
elseif isempty(maxIter)
    maxIter = Inf;
end

if nargin < 4
    s_interSteps = 25;
elseif isempty(s_interSteps)
    s_interSteps = 25;
end

if nargin < 3
    s_max = 1;
elseif isempty(s_max)
    s_max = 1;
end

if nargin < 2
    s_min = 1e-3;
elseif isempty(s_min)
    s_min = 1e-3;
end

nCoils = mrxisetup.getNumberOfCoils;

% variables:
% s ... coil segment locations
% p ... coil center locations
% c ... coil pattern locations in the origin, along z and with unit radius
% r ... radius
% R ... rotation matrix

unitCoilPatterns = mrxisetup.getCoilPatternWithUnitRadius();
VS = getVoxelSensorRelation(mrxisetup, false); % static voxel-sensor-relationship, 3 x nSensors x nVoxel matrix
r_init = mrxisetup.getCoilRadii();

%%%%%%%%%%%%% for testing
if debug
    if nCoils < 2
        error('Debugging only working for nCoils >= 2');
    end
    % get two coil indices from current setup
    ri = round(quantile(1:nCoils,[0.33 0.67]));
%     ri = [1 2];
    
    % set radius span to explore
    pointsperdim = 10;
    rd_span = [0.001 0.04];
    rd = rd_span(1):(rd_span(2)-rd_span(1))/pointsperdim:rd_span(2);
    
    % get initial radii
    rdb = mrxisetup.getCoilRadii();
    
    % initialize loop variables
    testcost = zeros(pointsperdim);
    testcond = testcost;
    yy = 0;
    for iy = rd
        yy = yy + 1;
        xx = 0;
        for ix = rd
            xx = xx + 1;
            rdb(ri(1)) = ix;
            rdb(ri(2)) = iy;
            [L_curr, L_notNormed, froL_single] = getLeadfieldMatrixWithNewCoilRadii(mrxisetup, VS, rdb, normL_flag);
            testcost(yy,xx) = algocost( L_curr );
            testcond(yy,xx) = cond( L_curr );
        end
    end


    % visualization
    figure(122); clf;
    imagesc(rd,rd,testcond);
    daspect([1 1 1]);
    shading interp;
    ci = colorbar;
    ci.Label.String = 'Spectral condition';
    hold on;
    scatter(r_init(ri(1)), r_init(ri(2)),'filled','k');
    %
    figure(123); clf;
    imagesc(rd,rd,testcost);
    daspect([1 1 1]);
    shading interp;
    ci = colorbar;
    ci.Label.String = 'Frobenius condition';
    hold on;
    scatter(r_init(ri(1)), r_init(ri(2)),'filled','k');
end
%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reset coil radii
r = r_init;
r(r<lowerRadiusLimit) = lowerRadiusLimit;
mrxisetup.setCoilRadii(r);

tic
% initialize loop parameters
iter = 0;
s_search = logspace(log10(s_min), log10(s_max), s_interSteps);
cost_temp = zeros(s_interSteps, 1);

% compute current leadfield matrix
[L_curr, L_notNormed, froL_single] = getLeadfieldMatrixWithNewCoilRadii(mrxisetup, VS, r, normL_flag);

cost = algocost( L_curr );
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
    % store previous mrxisetup
    mrxisetup_prev = mrxisetup.copySetup;
    
    % compute Frobenius norm of L
    froL = norm(L_curr,'fro');
    
    % compute Pseudo
    pL_curr = pinv(L_curr);
    fropL = norm(pL_curr,'fro');
    
    % compute derivative of Frobenius condition number
    r_grad = updateActivationGradient(mrxisetup, unitCoilPatterns, VS, L_curr, pL_curr, froL, fropL, ...
        L_notNormed, froL_single, normL_flag);
    
    % normalize gradient
    r_grad = r_grad./norm(r_grad);
    if debug
        r_grad((1:nCoils~=ri(1))&(1:nCoils~=ri(2))) = 0;
    end

   % determine stepsize with largest cost reduction within s_min and s_max
    % using s_interSteps with logarithmically equal steps
    for si = 1:s_interSteps
        rnew = r-s_search(si)*r_grad;
        rnew(rnew<lowerRadiusLimit) = lowerRadiusLimit;
        L_curr = getLeadfieldMatrixWithNewCoilRadii(mrxisetup, VS, rnew, normL_flag);
        cost_temp(si) = algocost( L_curr );
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
        mrxisetup = mrxisetup_prev;
        if verbose
            fprintf('Smallest stepsize %.3e increases cost. Terminating optimization.\n', s_search(1));
        end
        break;
    % else update stepsize and costs
    else
        s = s_search(si-1);
        r = r - s*r_grad;
        r(r<lowerRadiusLimit) = lowerRadiusLimit;
        mrxisetup.setCoilRadii(r);
        [L_curr, L_notNormed, froL_single] = getLeadfieldMatrixWithNewCoilRadii(mrxisetup, VS, r, normL_flag);
        iter = iter+1;
        cost_prev = cost;
        cost = cost_temp(si-1);
    end
    
    %%%%%%%%%%%%% for testing
    if debug
        figure(122)
        scatter(r(ri(1)), r(ri(2)),'r');
        figure(123)
        scatter(r(ri(1)), r(ri(2)),'r');
        pause(1e-3);
    end
    %%%%%%%%%%%%%

    cost_red = 100*(cost_prev-cost)/cost_prev; % cost reduction in percent
    % print cost if demanded
    if verbose
        fprintf('Current cost: %.5e, stepsize: %.3e, iteration: %3.0f, cost reduction: %3.4f%%, runtime: %.1f\n', cost, s, iter-1, cost_red,toc);
    end
%     % break if cost reduction is too small ( < 0.01%)
%     cost_thresh = 0.00000001;
%     if cost_red < cost_thresh
%         fprintf('Cost reduction below %g%%. Aborting.\n',cost_thresh);
%         break;
%     end
end

mrxisetup_optim = mrxisetup.copySetup;
cost_optim = cost;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print final cost if demanded
if verbose
    fprintf('Final cost: %.5e, runtime: %.1f\n', cost_optim, toc);
end

%%%%%%%%%%%%% for testing
if debug
    figure(122)
    scatter(r(ri(1)), r(ri(2)),'filled','g');
    figure(123)
    scatter(r(ri(1)), r(ri(2)),'filled','g');
    pause(1e-3);
end
%%%%%%%%%%%%%
end

function [L, L_notNormed, froL_single] = getLeadfieldMatrixWithNewCoilRadii(mrxisetup, VS, r, normL_flag)
    nCoils = mrxisetup.getNumberOfCoils();
    nSensors = mrxisetup.getNumberOfSensors();
    nVoxel = mrxisetup.getNumberOfVoxels();
    L = zeros(nCoils*nSensors, nVoxel);
    
    mrxisetup.setCoilRadii(r);
    VC = getVoxelCoilRelation(mrxisetup, false);
    if ~normL_flag
        for ci = 1:nCoils
            si = (nSensors*(ci-1)+1):nSensors*ci;
            L(si,:) = squeeze(sum(bsxfun(@times,VC(:,ci,:), VS),1));
        end
        L_notNormed = [];
        froL_single = [];
    else
        froL_single = zeros(nCoils,1);
        L_notNormed = L;
        for ci = 1:nCoils
            si = (nSensors*(ci-1)+1):nSensors*ci;
            L_notNormed(si,:) = squeeze(sum(bsxfun(@times,VC(:,ci,:), VS),1));
            froL_single(ci) = norm(L_notNormed(si,:), 'fro');
            L(si,:) = L_notNormed(si,:)./froL_single(ci);
        end
    end
end

function cost = algocost( L, tol )
    if nargin == 4
        L = L./normest(L);
        pL = pinv( L ,tol);
    else
        pL = pinv(L);
    end
    cost = norm(L,'fro') * norm(pL,'fro');
end


function r_grad = updateActivationGradient(mrxisetup, unitCoilPatterns, VS, L, pL, froL, fropL, ...
    L_notNormed, froL_single, normL_flag)
    % returns the derivative of the Frobenius condition of L w.r.t. the
    % coil radii 
    % OUTPUT: nCoils vector with gradient of cond(L,'fro') w.r.t. the coil
    % radii

    nCoils = mrxisetup.getNumberOfCoils;
    nVoxel = mrxisetup.getNumberOfVoxels;
    nSensors = mrxisetup.getNumberOfSensors;
    cp = mrxisetup.coilData.coilpatternassignment;
    
    p = permute(mrxisetup.ROIData.voxels,[1 3 2]);
    
    H_grad = zeros(3,nCoils,nVoxel);
    r_grad = zeros(nCoils,1);
    
    pLpLT = pL*pL';
    npL_div_nL = fropL/froL;
    nL_div_npL = 1/npL_div_nL;
    
    for ci = 1:nCoils
        % get static parameter R*c
        Rc1 = repmat(mrxisetup.coilData.rotMatrices(:,:,ci)*unitCoilPatterns{cp(ci)}(:,1:end), [1 1 nVoxel]);
        Rc2 = circshift(Rc1, -1, 2);
        % get segment vectors s of current coil
        s1 = bsxfun(@minus, mrxisetup.coilData.segments.(['c' num2str(ci)])(:,1:end), p);
        s2 = circshift(s1, -1, 2);
        % compute norms of s squared
        ns_sq = sum(s1.*s1, 1);
        % compute norms of s
        ns1 = sqrt(ns_sq);
        ns2 = circshift(ns1, -1, 2);
        % compute product of norms ||s1||*||s2||
        ns1ns2 = ns1.*ns2;
        % compute scalar products <s, Rc>
        SP_s_Rc = dot(s1, Rc1, 1);
        % compute scalar products <s, Rc> divided by ||s||
        SP_s_Rc_div_ns1 = SP_s_Rc./ns1;
        SP_s_Rc_div_ns2 = circshift(SP_s_Rc_div_ns1, -1, 2);
        % compute cross products s1 x s2
        CP_s1_s2 = cross(s1, s2, 1);
        % compute scalar products <s1, s2>
        SP_s1_s2 = dot(s1, s2, 1);
        
        %%% compute terms of the derivative of H step by step
        % numerator
        Num =   bsxfun(@times, (ns1 + ns2), CP_s1_s2);
        % derivative of numerator
        dNum =  bsxfun(@times, (SP_s_Rc_div_ns1 + SP_s_Rc_div_ns2), CP_s1_s2) +...
                bsxfun(@times, (ns1 + ns2), cross(s1, Rc2, 1) + cross(Rc1, s2, 1));
        % denominator
        Den =   ns1ns2 .* (ns1ns2 + SP_s1_s2);
        % derivative of denominator
        dDen =  2*(ns_sq.*circshift(SP_s_Rc, -1, 2) + circshift(ns_sq, -1, 2).*SP_s_Rc) +...
                (SP_s_Rc_div_ns1.*ns2 + SP_s_Rc_div_ns2.*ns1).*SP_s1_s2 +...
                ns1ns2.*(dot(Rc1, s2, 1)+dot(Rc2, s1, 1));

        %%% compute gradient of H
        H_grad(:,ci,:) =...
            sum(bsxfun(@rdivide,...
            (bsxfun(@times, dNum(:,1:end-1,:), Den(:,1:end-1,:)) - bsxfun(@times, Num(:,1:end-1,:), dDen(:,1:end-1,:))),...
            4*pi*Den(:,1:end-1,:).^2),2); 
        
        %%% compute L_grad = L(H_grad)
        L_grad = squeeze(sum(bsxfun(@times,H_grad(:,ci,:), VS),1));
        
        % compute gradient w.r.t. normalized leadfield matrix if
        % demanded
        ii = nSensors*(ci-1)+1:nSensors*ci;
        if normL_flag
            L_grad =    L_grad./froL_single(ci) - ...
                        sum(sum(L_grad.*L_notNormed(ii,:)))/froL_single(ci)^3 * L_notNormed(ii,:);
        end
        
        %%% compute r_grad = d(cond(L,'fro'))/dr
        ii = nSensors*(ci-1)+1:nSensors*ci;
        r_grad(ci) = sum(sum(L(ii,:).*L_grad))*npL_div_nL - nL_div_npL*sum(sum(pLpLT.*(pL(:,ii)*L_grad)));
    end
end
