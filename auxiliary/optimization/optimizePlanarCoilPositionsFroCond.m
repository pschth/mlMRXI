function [ mrxisetup_optim, cost_optim] = optimizePlanarCoilPositionsFroCond( ...
    mrxisetup_init, constraintAlloc, s_min, s_max, s_interSteps, maxIter, normL_flag, verbose, debug )
% optimizes the coil currents in MRXI by minimizing the Frobenius condition number of
% the MRXI Leadfield Matrix

% INPUT
% mrxisetup ... standard MRXI SETUP object
% constraintAlloc ... nCoils vector; contains for each coil the constraint
% index the respective coil is assigned to, Default: ones(nCoils,1)
% s_min (optional) ... minimum stepsize, Default: 1e-3
% s_max (optional) ... maximum stepsize, Default: 1
% s_interSteps (optional) ... number of logarithmically equidistant steps
% between s_min and s_max, Default: 25
% maxIter (optional) ... maximum number of iterations (per activation if
% normL_flag (optional) ... boolean, if true, the individual system
% matrices of each coil are normalized such that Lc = Lc./norm(Lc,'fro').
% Default: false
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
    normL_flag = false;
elseif isempty(normL_flag)
    normL_flag = false;
end

if nargin < 6
    maxIter = Inf;
elseif isempty(maxIter)
    maxIter = Inf;
end

if nargin < 5
    s_interSteps = 25;
elseif isempty(s_interSteps)
    s_interSteps = 25;
end

if nargin < 4
    s_max = 1;
elseif isempty(s_max)
    s_max = 1;
end

if nargin < 3
    s_min = 1e-3;
elseif isempty(s_min)
    s_min = 1e-3;
end

if nargin < 2
    constraintAlloc = ones(mrxisetup.getNumberOfCoils,1);
elseif isempty(s_min)
    constraintAlloc = ones(mrxisetup.getNumberOfCoils,1);
end

% check if planar constraint(s) was/were defined, otherwise abort
if ~isfield(mrxisetup.coilData, 'constraints')
    error('MRXI setup object has no defined constraints.');
else
    for cp = 1:length(mrxisetup.coilData.constraints)
        if ~strcmp(mrxisetup.coilData.constraints{cp}.type, 'planar')
            error('Only planar constraints allowed.');
        end
    end
end


VS = getVoxelSensorRelation(mrxisetup, false); % static voxel-sensor-relationship, 3 x nSensors x nVoxel matrix

%%%%%%%%%%%%%% for testing
% works only for 1 coil and 1 planar constraint
if debug
    % disable debugging if conditions aren't met
    if mrxisetup.getNumberOfCoils ~= 1 || length(mrxisetup.coilData.constraints) ~= 1
       warning('Debugging works only for 1 coil and 1 planar coil constraint. Disabling debugging mode.');
       debug=false;
    elseif ~strcmp(mrxisetup.coilData.constraints{1}.type, 'planar')
       warning('Debugging works only for 1 coil and 1 planar coil constraint. Disabling debugging mode.');
       debug=false;
    end
end
if debug
    % copy setup
    mrxisetup_debug = mrxisetup.copySetup;
    % resolution per dimension of cost computation on plane
    resPerDim = 20;
    s = linspace(0,1,resPerDim);
    interpPoints = zeros(3, resPerDim^2);
    a1 = mrxisetup_debug.coilData.constraints{1}.a1;
    a2 = mrxisetup_debug.coilData.constraints{1}.a2;

    % generate grid points for cost computation
    for r = 1:resPerDim
       currIdx = (r-1)*resPerDim+1:r*resPerDim;
       interpPoints(:,currIdx) = bsxfun(@plus, s(r)*a2, bsxfun(@times, s, a1));
    end
    interpPoints = bsxfun(@plus, interpPoints, mrxisetup_debug.coilData.constraints{1}.a0);

    % compute costs of grid points
    cost_debug = zeros(1,resPerDim^2);
    for r = 1:resPerDim^2
       cost_debug(r) = algocost( getLeadfieldMatrixWithNewCoilPositions(mrxisetup_debug, VS, interpPoints(:,r)) );
    end

    % visualize plane
    X = reshape(interpPoints(1,:), resPerDim, resPerDim);
    Y = reshape(interpPoints(2,:), resPerDim, resPerDim);
    Z = reshape(interpPoints(3,:), resPerDim, resPerDim);
    C = reshape(cost_debug, resPerDim, resPerDim);
    figure(121); clf;
    surf(X,Y,Z,C);
    shading interp;
    hold on

    % visualize initial coil position
    scatter3(mrxisetup.coilData.centers(1,1),...
            mrxisetup.coilData.centers(2,1),...
            mrxisetup.coilData.centers(3,1),'filled','k');
    
    % visualize logarithmic cost function for better supervision
    figure(122); clf;
    surf(X,Y,Z,log10(C));
    shading interp;
    hold on
    scatter3(mrxisetup.coilData.centers(1,1),...
            mrxisetup.coilData.centers(2,1),...
            mrxisetup.coilData.centers(3,1),'filled','k');
    drawnow;
end
%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
% initialize loop parameters
iter = 0;
s_search = logspace(log10(s_min), log10(s_max), s_interSteps);
current_displacement_vectors = zeros(3,mrxisetup.getNumberOfCoils);
cost_temp = zeros(s_interSteps, 1);

% get current coil centers and compute current leadfield matrix
centers = mrxisetup.coilData.centers;
[L_curr, L_notNormed, froL_single] = getLeadfieldMatrixWithNewCoilPositions(mrxisetup, VS, centers, normL_flag);

% check allocated constraints
usedConstraints = unique(constraintAlloc);
if max(usedConstraints) > length(mrxisetup.coilData.constraints)
    error('Specified constraintAlloc index(es) exceed number of constraints.');
else
    % else, precompute matrices for gradient calculation
    nCoilPatterns = length(mrxisetup.coilData.coilpattern);
    nSegments = zeros(nCoilPatterns,1);
    for cp = 1:nCoilPatterns
        nSegments(cp) = size(mrxisetup.coilData.coilpattern{cp},2);
    end
    nVoxels = mrxisetup.getNumberOfVoxels;
    nCoils = mrxisetup.getNumberOfCoils;
    anMatrix = cell(2,length(usedConstraints),nCoilPatterns);
    for constraint = usedConstraints(:)'
        for cp = 1:nCoilPatterns
            anMatrix{1,constraint,cp} = repmat(mrxisetup.coilData.constraints{constraint}.a1n,...
                                            1, nSegments(mrxisetup.coilData.coilpatternassignment(cp)), nVoxels);
            anMatrix{2,constraint,cp} = repmat(mrxisetup.coilData.constraints{constraint}.a2n,...
                                            1, nSegments(mrxisetup.coilData.coilpatternassignment(cp)), nVoxels);
        end
    end
end

% get coil extensions for constraint check
coilExtensions = mrxisetup.getPlanarCoilExtensions;
% get constraint extensions for constraint check
nConstraints = length(usedConstraints);
constraintExtensions = zeros(2,nConstraints);
for constraint = 1:nConstraints
    currConstraint = usedConstraints(constraint);
    constraintExtensions(:,constraint) = [  norm(mrxisetup.coilData.constraints{currConstraint}.a1); ...
                                            norm(mrxisetup.coilData.constraints{currConstraint}.a2)];
end

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
    ccenter_grad = updateActivationGradient(...
        mrxisetup, constraintAlloc, anMatrix, VS, ...
        L_curr, pL_curr, froL, fropL, froL_single, L_notNormed);
    
    % normalize gradient
    ccenter_grad = ccenter_grad./norm(ccenter_grad,'fro');

   % determine stepsize with largest cost reduction within s_min and s_max
    % using s_interSteps with logarithmically equidistant steps
    for si = 1:s_interSteps
        for ci = 1:nCoils
            constraint = constraintAlloc(ci);
            current_displacement_vectors(:,ci) = s_search(si) * (...
                mrxisetup.coilData.constraints{constraint}.a1n * ccenter_grad(1,ci) + ...
                mrxisetup.coilData.constraints{constraint}.a2n * ccenter_grad(2,ci));
        end
        
        centers_temp = centers - current_displacement_vectors;
        
        % check if coil is still inside borders of constraint - if not,
        % move to border
        for ci = 1:nCoils
            % get extension from current constraint
            constraint = constraintAlloc(ci);
            currConstraintExtensions = constraintExtensions(:, usedConstraints == constraint);
            
            % compute location on constraint and current coil extensions
            locOnConstraint = centers_temp(:,ci) - mrxisetup.coilData.constraints{constraint}.a0;
            currCoilExtensions = coilExtensions(:,ci)./2;
            
            % compute projections of maximum coil extensions on constraint
            % borders
            projOna1n = locOnConstraint'*mrxisetup.coilData.constraints{constraint}.a1n;
            projOna2n = locOnConstraint'*mrxisetup.coilData.constraints{constraint}.a2n;
            projMaxOna1n = projOna1n + currCoilExtensions(1);
            projMinOna1n = projOna1n - currCoilExtensions(1);
            projMaxOna2n = projOna2n + currCoilExtensions(2);
            projMinOna2n = projOna2n - currCoilExtensions(2);
            
            %%%% check if maximum/minimum coil extensions are inside border
            %%%% - else move
            % check a1 max and a1 min
            if projMaxOna1n > currConstraintExtensions(1)
                dist = projMaxOna1n - currConstraintExtensions(1);
                centers_temp(:,ci) = centers_temp(:,ci) - dist * mrxisetup.coilData.constraints{constraint}.a1n;
            elseif projMinOna1n < 0
                dist = 0 - projMinOna1n;
                centers_temp(:,ci) = centers_temp(:,ci) + dist * mrxisetup.coilData.constraints{constraint}.a1n;
            end
            % check a2 max and a2 min
            if projMaxOna2n > currConstraintExtensions(2)
                dist = projMaxOna2n - currConstraintExtensions(2);
                centers_temp(:,ci) = centers_temp(:,ci) - dist * mrxisetup.coilData.constraints{constraint}.a2n;
            elseif projMinOna2n < 0
                dist = 0 - projMinOna2n;
                centers_temp(:,ci) = centers_temp(:,ci) + dist * mrxisetup.coilData.constraints{constraint}.a2n;
            end
        end
        
        L_curr = getLeadfieldMatrixWithNewCoilPositions(mrxisetup, VS, centers_temp, normL_flag);
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
        for ci = 1:nCoils
            constraint = constraintAlloc(ci);
            current_displacement_vectors(:,ci) = s * (...
                mrxisetup.coilData.constraints{constraint}.a1n * ccenter_grad(1,ci) + ...
                mrxisetup.coilData.constraints{constraint}.a2n * ccenter_grad(2,ci));
        end
        % compute new centers
        centers = centers - current_displacement_vectors;

        % check if coil is still inside borders of constraint - if not,
        % move to border
        for ci = 1:nCoils
            % get extension from current constraint
            constraint = constraintAlloc(ci);
            currConstraintExtensions = constraintExtensions(:, usedConstraints == constraint);
            
            % compute location on constraint and current coil extensions
            locOnConstraint = centers(:,ci) - mrxisetup.coilData.constraints{constraint}.a0;
            currCoilExtensions = coilExtensions(:,ci)./2;
            
            % compute projections of maximum coil extensions on constraint
            % borders
            projOna1n = locOnConstraint'*mrxisetup.coilData.constraints{constraint}.a1n;
            projOna2n = locOnConstraint'*mrxisetup.coilData.constraints{constraint}.a2n;
            projMaxOna1n = projOna1n + currCoilExtensions(1);
            projMinOna1n = projOna1n - currCoilExtensions(1);
            projMaxOna2n = projOna2n + currCoilExtensions(2);
            projMinOna2n = projOna2n - currCoilExtensions(2);
            
            %%%% check if maximum/minimum coil extensions are inside border
            %%%% - else move
            % check a1 max and a1 min
            if projMaxOna1n > currConstraintExtensions(1)
                dist = projMaxOna1n - currConstraintExtensions(1);
                centers(:,ci) = centers(:,ci) - dist * mrxisetup.coilData.constraints{constraint}.a1n;
            elseif projMinOna1n < 0
                dist = 0 - projMinOna1n;
                centers(:,ci) = centers(:,ci) + dist * mrxisetup.coilData.constraints{constraint}.a1n;
            end
            % check a2 max and a2 min
            if projMaxOna2n > currConstraintExtensions(2)
                dist = projMaxOna2n - currConstraintExtensions(2);
                centers(:,ci) = centers(:,ci) - dist * mrxisetup.coilData.constraints{constraint}.a2n;
            elseif projMinOna2n < 0
                dist = 0 - projMinOna2n;
                centers(:,ci) = centers(:,ci) + dist * mrxisetup.coilData.constraints{constraint}.a2n;
            end
        end
        
        mrxisetup.setCoilCenters(centers);
        [L_curr, L_notNormed, froL_single] = getLeadfieldMatrixWithNewCoilPositions(mrxisetup, VS, centers, normL_flag);
        iter = iter+1;
        cost_prev = cost;
        cost = cost_temp(si-1);
    end
    
    %%%%%%%%%%%%%% for testing
    if debug
        figure(121);
        % visualize optimizitation steps of coil positioning
        scatter3(mrxisetup.coilData.centers(1,1),...
                mrxisetup.coilData.centers(2,1),...
                mrxisetup.coilData.centers(3,1),'r');
        drawnow;
    end
    %%%%%%%%%%%%%%

    cost_red = 100*(cost_prev-cost)/cost_prev; % cost reduction in percent
    % print cost if demanded
    if verbose
        fprintf('Current cost: %.5e, stepsize: %.3e, iteration: %3.0f, cost reduction: %3.4f%%, runtime: %.1f\n', cost, s, iter-1, cost_red,toc);
    end
%     % break if cost reduction is too small ( < 0.01%)
%     cost_thresh = 0.01;
%     if cost_red < cost_thresh
%         fprintf('Cost reduction below %g%%. Aborting.\n',cost_thresh);
%         break;
%     end
end

%%%%%%%%%%%%%% for testing
if debug
    % visualize optimized coil position
    figure(121);
    scatter3(mrxisetup.coilData.centers(1,1),...
            mrxisetup.coilData.centers(2,1),...
            mrxisetup.coilData.centers(3,1),'filled','g');
    figure(122);
    scatter3(mrxisetup.coilData.centers(1,1),...
            mrxisetup.coilData.centers(2,1),...
            mrxisetup.coilData.centers(3,1),'filled','g');
    drawnow;
end
%%%%%%%%%%%%%%

mrxisetup_optim = mrxisetup.copySetup;
cost_optim = cost;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print final cost if demanded
if verbose
    fprintf('Final cost: %.5e, runtime: %.1f\n', cost_optim, toc);
end

end

function [L, L_notNormed, froL_single] = getLeadfieldMatrixWithNewCoilPositions(mrxisetup, VS, centers, normL_flag)
    nCoils = mrxisetup.getNumberOfCoils();
    nSensors = mrxisetup.getNumberOfSensors();
    nVoxel = mrxisetup.getNumberOfVoxels();
    L = zeros(nCoils*nSensors, nVoxel);
    if normL_flag
        L_notNormed = L;
        froL_single = zeros(nCoils,1);
    else
        L_notNormed = [];
        froL_single = [];
    end
    
    mrxisetup.setCoils(centers, mrxisetup.coilData.orientations, ...
        mrxisetup.coilData.coilpattern, ...
        mrxisetup.coilData.coilpatternassignment);
    mrxisetup.setCoilRadii(mrxisetup.getCoilRadii);
    
    VC = getVoxelCoilRelation(mrxisetup, false);
    for ci = 1:nCoils
        si = (nSensors*(ci-1)+1):nSensors*ci;
        if normL_flag
            L_notNormed(si,:) = squeeze(sum(bsxfun(@times,VC(:,ci,:), VS),1));
            froL_single(ci) = norm(L_notNormed(si,:),'fro');
            L(si,:) = L_notNormed(si,:)./froL_single(ci);
        else
            L(si,:) = squeeze(sum(bsxfun(@times,VC(:,ci,:), VS),1));
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


function ccenter_grad = updateActivationGradient(mrxisetup, constraintAlloc, anMatrix, VS, L, pL, froL, fropL, froL_single, L_notNormed)
    % returns the derivative of the Frobenius condition of L w.r.t. the
    % coil positions on a planar constraint 
    % OUTPUT: nCoils vector with gradient of cond(L,'fro') w.r.t. the coil
    % positions

    nCoils = mrxisetup.getNumberOfCoils;
    nSensors = mrxisetup.getNumberOfSensors;
    normL_flag = ~isempty(froL_single);
    
    p = permute(mrxisetup.ROIData.voxels,[1 3 2]);
    
    ccenter_grad = zeros(2,nCoils);
    
    pLpLT = pL*pL';
    npL_div_nL = fropL/froL;
    nL_div_npL = 1/npL_div_nL;
    
    for ci = 1:nCoils
        % get current plane vectors for the right constraint
        constraint = constraintAlloc(ci);
        cp = mrxisetup.coilData.coilpatternassignment(ci);
        an(:,1) = mrxisetup.coilData.constraints{constraint}.a1n; % = ds1
        an(:,2) = mrxisetup.coilData.constraints{constraint}.a2n; % = ds2
        % get segment vectors s of current coil
        s1 = bsxfun(@minus, mrxisetup.coilData.segments.(['c' num2str(ci)])(:,1:end), p);
        s2 = circshift(s1, -1, 2);
        % compute norms of s squared
        ns_sq1 = sum(s1.*s1, 1);
        ns_sq2 = circshift(ns_sq1, -1, 2);
        % compute norms of s
        ns1 = sqrt(ns_sq1);
        ns2 = circshift(ns1, -1, 2);
        % compute sum of norms
        ns1_p_ns2 = ns1 + ns2;
        % compute difference of vectors
        s2_m_s1 = s2-s1;
        % compute product of norms ||s1||*||s2||
        ns1ns2 = ns1.*ns2;
        % compute s/ns and then s1/ns1 + s2/ns2
        s1_d_ns1 = bsxfun(@rdivide, s1, ns1);
        s2_d_ns2 = circshift(s1_d_ns1, -1, 2);
        s1dns1_p_s2dns2 = s1_d_ns1 + s2_d_ns2;
        % compute cross products s1 x s2
        CP_s1_s2 = cross(s1, s2, 1);
        % compute scalar products <s1, s2>
        SP_s1_s2 = dot(s1, s2, 1);

        %%% compute terms of the derivative of H step by step
        % numerator
        Num =   bsxfun(@times, ns1_p_ns2, CP_s1_s2);
        % denominator
        Den =   ns1ns2 .* (ns1ns2 + SP_s1_s2);
        % derivative of denominator without scalar product with ds
        dDen_wo_ds =    bsxfun(@times, 2*ns_sq2 + ns2./ns1.*SP_s1_s2 + ns1ns2, s1) + ...
                        bsxfun(@times, 2*ns_sq1 + ns1./ns2.*SP_s1_s2 + ns1ns2, s2);

        for m = 1:2
            % derivative of numerator
            dNum =  bsxfun(@times, sum(bsxfun(@times, an(:,m), s1dns1_p_s2dns2), 1), CP_s1_s2) +...
                    bsxfun(@times, ns1_p_ns2, cross( anMatrix{m,constraint,cp}, s2_m_s1, 1) );
            % derivative of denominator
            dDen =  sum(bsxfun(@times, an(:,m), dDen_wo_ds),1);

            %%% compute gradient of H
            H_grad =...
                sum(bsxfun(@rdivide,...
                (bsxfun(@times, dNum(:,1:end-1,:), Den(:,1:end-1,:)) - bsxfun(@times, Num(:,1:end-1,:), dDen(:,1:end-1,:))),...
                4*pi*Den(:,1:end-1,:).^2),2); 

            %%% compute L_grad = L(H_grad)
            L_grad = squeeze(sum(bsxfun(@times,H_grad, VS),1));
            % compute gradient w.r.t. normalized leadfield matrix if
            % demanded
            ii = nSensors*(ci-1)+1:nSensors*ci;
            if normL_flag
                L_grad =    L_grad./froL_single(ci) - ...
                            sum(sum(L_grad.*L_notNormed(ii,:)))/froL_single(ci)^3 * L_notNormed(ii,:);
            end
            %%% compute ccenter_grad = d(cond(L,'fro'))/dr
            ii = nSensors*(ci-1)+1:nSensors*ci;
            ccenter_grad(m,ci) = sum(sum(L(ii,:).*L_grad))*npL_div_nL - nL_div_npL*sum(sum(pLpLT.*(pL(:,ii)*L_grad)));
        end
    end
end
