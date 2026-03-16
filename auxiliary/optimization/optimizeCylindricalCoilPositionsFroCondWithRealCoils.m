function [ mrxisetup_optim, cost_optim] = ...
    optimizeCylindricalCoilPositionsFroCondWithRealCoils( ...
    mrxisetup_init, constraintAlloc, s_min, s_max, s_interSteps, maxIter, normL_flag, verbose, debug )
% optimizes the coil positions in MRXI by minimizing the Frobenius
% condition number of the MRXI Leadfield Matrix along the cylindrical
% constraints defined

% INPUT
% mrxisetup_init ... standard MRXI SETUP object from where the optimization
% starts
% constraintAlloc ... nCoils vector; contains for each coil the constraint
% index the respective coil is assigned to, Default: ones(nCoils,1)
% s_min (optional) ... minimum stepsize, Default: 1e-3
% s_max (optional) ... maximum stepsize, Default: 1
% s_interSteps (optional) ... number of logarithmically equidistant steps
% between s_min and s_max, Default: 25
% maxIter (optional) ... maximum number of iterations, Default: 100
% normL_flag (optional) ... boolean, if true, the individual system
% matrices of each coil are normalized such that Lc = Lc./norm(Lc,'fro').
% Default: false
% verbose (optional) ... boolean, if false, no output is given from the function,
% Default: true
% debug (optional) ... boolean, if true, the optimization steps are
% visualized. Works only for 1 coil, otherwise it gets deactivated.
% Default: false

% OUTPUT
% mrxisetup_optim ... standard MRXI SETUP object; optimized coil setup
% cost_optim ... value of the cost function using optimized coil currents

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
    maxIter = 100;
elseif isempty(maxIter)
    maxIter = 100;
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

% check if cylindrical constraint(s) was/were defined, otherwise abort
if ~isfield(mrxisetup.coilData, 'constraints')
    error('MRXI setup object has no defined constraints.');
else
    for cp = 1:length(mrxisetup.coilData.constraints)
        if ~strcmp(mrxisetup.coilData.constraints{cp}.type, 'cylindrical')
            error('Only cylindrical constraints allowed.');
        end
    end
end


VS = getVoxelSensorRelation(mrxisetup, false); % static voxel-sensor-relationship, 3 x nSensors x nVoxel matrix
nCoils = mrxisetup.getNumberOfCoils;
% get maximum cylinder height limitation for each coil
hmax = zeros(nCoils,1);
for ci = 1:nCoils
    constIdx = constraintAlloc(ci);
    hmax(ci) = norm(mrxisetup.coilData.constraints{constIdx}.topCoord - ...
                mrxisetup.coilData.constraints{constIdx}.bottomCoord);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%% for debugging
% ATTENTION: works only for 1 coil!
if debug
    if mrxisetup.getNumberOfCoils > 1
        warning('debugging does not work for more than 1 coil. Disabling debug mode.')
        debug = false;
    else
        % copy setup
        mrxisetup_debug = mrxisetup.copySetup;
        % resolution per dimension of cost computation on plane
        resPerDim = 60;
        s = linspace(0,1,resPerDim);
        interpCenters = zeros(3, resPerDim^2);
        interpOrientations = interpCenters;
        a1 = mrxisetup_debug.coilData.constraints{constraintAlloc(1)}.bottomCoord;
        a2 = mrxisetup_debug.coilData.constraints{constraintAlloc(1)}.topCoord;
        r = mrxisetup_debug.coilData.constraints{constraintAlloc(1)}.radius;
        cAxis = a2-a1;

        % get rotation matrix to constraint axis
        [~, rotConst] = rotatepoints([0 0 1]', cAxis);
    %     rotMatrices = getRotationMatricesFromInitPosition(angles, anglesInit, constraintData)
        % generate grid points for cost computation
        i = 0;
        for r1 = 1:resPerDim
            for r2 = 1:resPerDim
               i = i+1;
               curr_height = s(r1) * norm(cAxis);
               curr_angle = s(r2) * 2*pi;
               unrotatedPos = [r*cos(curr_angle); r*sin(curr_angle); curr_height];
               interpCenters(:,i) = a1 + rotConst*unrotatedPos;
               interpOrientations(:,i) = normcols(-rotConst*[unrotatedPos(1:2);0]);
            end
        end

        % compute costs of grid points
        cost_debug = zeros(1,resPerDim^2);
        for ri = 1:resPerDim^2
           cost_debug(ri) = algocost( getLeadfieldMatrixWithNewCoilPositions(...
               mrxisetup_debug, VS, interpCenters(:,ri), interpOrientations(:,ri), normL_flag ));
        end

        % visualize plane
        X = reshape(interpCenters(1,:), resPerDim, resPerDim);
        Y = reshape(interpCenters(2,:), resPerDim, resPerDim);
        Z = reshape(interpCenters(3,:), resPerDim, resPerDim);
        C = reshape(cost_debug, resPerDim, resPerDim);
        figure(121); clf;
        surf(X,Y,Z,C);
        shading interp;
        hold on

        % visualize initial coil position
        quiverScale = min(diff(mrxisetup.getTransformedROIBoundary,[],2))/2;
        scatter3(mrxisetup.coilData.centers(1,1),...
                mrxisetup.coilData.centers(2,1),...
                mrxisetup.coilData.centers(3,1),'filled','r');
        quiver3(mrxisetup.coilData.centers(1,1),...
                mrxisetup.coilData.centers(2,1),...
                mrxisetup.coilData.centers(3,1),...
                mrxisetup.coilData.orientations(1,1)*quiverScale,...
                mrxisetup.coilData.orientations(2,1)*quiverScale,...
                mrxisetup.coilData.orientations(3,1)*quiverScale,'r');
        daspect([1 1 1]);

        % visualize logarithmic cost function for better supervision
        figure(122); clf;
        surf(X,Y,Z,log10(C));
        shading interp;
        hold on
        scatter3(mrxisetup.coilData.centers(1,1),...
                mrxisetup.coilData.centers(2,1),...
                mrxisetup.coilData.centers(3,1),'filled','r');
        quiver3(mrxisetup.coilData.centers(1,1),...
            mrxisetup.coilData.centers(2,1),...
            mrxisetup.coilData.centers(3,1),...
            mrxisetup.coilData.orientations(1,1)*quiverScale,...
            mrxisetup.coilData.orientations(2,1)*quiverScale,...
            mrxisetup.coilData.orientations(3,1)*quiverScale,'r');
        drawnow;
        daspect([1 1 1]);
    end
end
%% %%%%%%%%%%%%%%

tic
% initialize loop parameters
iter = 0;
s_search = logspace(log10(s_min), log10(s_max), s_interSteps);
cost_temp = zeros(s_interSteps, 1);

% get current coil centers and compute current leadfield matrix
centers = mrxisetup.coilData.centers;
orientations = mrxisetup.coilData.orientations;
[L_curr, L_notNormed, froL_single] = getLeadfieldMatrixWithNewCoilPositions(mrxisetup, VS, centers, orientations, normL_flag);

% check allocated constraints
constraintData.constraintAlloc =  constraintAlloc;
constraintData.usedConstraints = unique(constraintAlloc);
nUsedConstraints = length(constraintData.usedConstraints);

if max(constraintData.usedConstraints) > length(mrxisetup.coilData.constraints)
    error('Specified constraintAlloc index(es) exceed number of constraints.');
else
    % else, precompute matrices for gradient calculation
    nCoilPatterns = length(mrxisetup.coilData.coilpattern);
    nSegments = zeros(nCoilPatterns,1);
    for cp = 1:nCoilPatterns
        nSegments(cp) = size(mrxisetup.coilData.coilpattern{cp},2);
    end
    % compute rotation matrices, central axes as well as normal axes where
    % angle = 0 for every constraint
    constraintData.rotMatrixConst = zeros(3,3,nUsedConstraints);
    constraintData.axesConst = zeros(3,nUsedConstraints);
    constraintData.angleRefAxesConst = constraintData.axesConst;
    constraintData.angleDirAxesConst = constraintData.axesConst;
    for ci = 1:nUsedConstraints
        % get index of used constraint
        cIdx = constraintData.usedConstraints(ci);
        % compute central axes
        constraintData.axesConst(:,ci) = mrxisetup.coilData.constraints{cIdx}.topCoord - ...
                                         mrxisetup.coilData.constraints{cIdx}.bottomCoord;
        % compute rotation matrices
        [~, constraintData.rotMatrixConst(:,:,ci)] = rotatepoints([0 0 1]', constraintData.axesConst(:,ci));
        % define reference axis for angle offset
        constraintData.angleRefAxesConst(:,ci) = constraintData.rotMatrixConst(:,:,ci) * [1 0 0]';
        % define directional axis for angle offset
        constraintData.angleDirAxesConst(:,ci) = constraintData.rotMatrixConst(:,:,ci) * [0 1 0]';
    end
end


[anglesOnConst_INIT, heightsCurr] = getCylindricalCoords(mrxisetup, constraintData);
anglesCurr = anglesOnConst_INIT;

cost = algocost( L_curr );

%%%%%%%%%%%%%%%%%%% debugging:

if debug
    [H,A] = meshgrid(linspace(0,norm(cAxis),resPerDim), linspace(0,2*pi,resPerDim));
    figure(123); clf;
    surf(H,A,log10(C));
    shading interp;
    hold on
    scatter3(heightsCurr, anglesCurr, log10(cost),'filled','r');
end

%%%%%%%%%%%%%%%%%%%

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
    coilGrad = updatePositionGradient(...
        mrxisetup, mrxisetup_init, constraintData, anglesOnConst_INIT, anglesCurr, ...
        VS, L_curr, pL_curr, froL, fropL, froL_single, L_notNormed);
    
    % normalize gradient
    coilGrad = coilGrad./norm(coilGrad,'fro')*nCoils;
%     coilGrad = pinv(coilGrad'); % = gauss-newton steps instead of gradient descent
    %%%%%%%%%%%%%%%%%%% debugging:

    if debug
        figure(123);
        quiver3(heightsCurr, anglesCurr, log10(cost), ...
                -coilGrad(1)*quiverScale, -coilGrad(2)*quiverScale, 0, 'k');
        drawnow;
    end

    %%%%%%%%%%%%%%%%%%%

   % determine stepsize with largest cost reduction within s_min and s_max
    % using s_interSteps with logarithmically equidistant steps
    cost_temp(:) = 0;
    for si = 1:s_interSteps
        
        anglesNext = anglesCurr - s_search(si)*coilGrad(2,:)';
        heightNext = heightsCurr- s_search(si)*coilGrad(1,:)';
        % check if border limitations are surpassed
        for ci = 1:nCoils
            if heightNext(ci) < 0
                heightNext(ci) = 0;
            elseif heightNext(ci) > hmax(ci)
                heightNext(ci) = hmax(ci);
            end
        end
        [centers, orientations] = ...
            getEuclideanCoordinates(mrxisetup_init, anglesOnConst_INIT, ...
                                    anglesNext, ...
                                    heightNext, ...
                                    constraintData);
        
        L_curr = getLeadfieldMatrixWithNewCoilPositions(mrxisetup, VS, centers, orientations, normL_flag);
        cost_temp(si) = algocost( L_curr );
        lastFlag = false;
        if si > 1
            if cost_temp(si) > cost_temp(si-1)
                break;
            elseif si == s_interSteps
                lastFlag = true;
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
        if lastFlag
            s = s_search(si);
        else
            s = s_search(si-1);
        end
        anglesCurr = mod(anglesCurr - s*coilGrad(2,:)',2*pi);
        heightsCurr = heightsCurr -s*coilGrad(1,:)';
        % check if border limitations are surpassed
        for ci = 1:nCoils
            if heightsCurr(ci) < 0
                heightsCurr(ci) = 0;
            elseif heightsCurr(ci) > hmax(ci)
                heightsCurr(ci) = hmax(ci);
            end
        end
        [centers, orientations] = getEuclideanCoordinates(mrxisetup_init, anglesOnConst_INIT, anglesCurr, heightsCurr, constraintData);
        R = mrxisetup.getCoilRadii;
        mrxisetup.setCoils(centers, orientations, ...
            mrxisetup.coilData.coilpattern, ...
            mrxisetup.coilData.coilpatternassignment);
        mrxisetup.setCoilRadii(R);
        [L_curr, L_notNormed, froL_single] = getLeadfieldMatrixWithNewCoilPositions(mrxisetup, VS, centers, orientations, normL_flag);
        iter = iter+1;
        cost_prev = cost;
        cost = cost_temp(si-1);
        
        %%%%%%%%%%%%% debugging
        if debug
            figure(121);
            scatter3(mrxisetup.coilData.centers(1,1),...
                    mrxisetup.coilData.centers(2,1),...
                    mrxisetup.coilData.centers(3,1),'k*');
            quiver3(mrxisetup.coilData.centers(1,1),...
                mrxisetup.coilData.centers(2,1),...
                mrxisetup.coilData.centers(3,1),...
                mrxisetup.coilData.orientations(1,1)*quiverScale,...
                mrxisetup.coilData.orientations(2,1)*quiverScale,...
                mrxisetup.coilData.orientations(3,1)*quiverScale,'k');
            
            figure(122);
            scatter3(mrxisetup.coilData.centers(1,1),...
                    mrxisetup.coilData.centers(2,1),...
                    mrxisetup.coilData.centers(3,1),'k*');
            quiver3(mrxisetup.coilData.centers(1,1),...
                mrxisetup.coilData.centers(2,1),...
                mrxisetup.coilData.centers(3,1),...
                mrxisetup.coilData.orientations(1,1)*quiverScale,...
                mrxisetup.coilData.orientations(2,1)*quiverScale,...
                mrxisetup.coilData.orientations(3,1)*quiverScale,'k');
            
            figure(123);
            scatter3(heightsCurr, anglesCurr, log10(cost), 'k');
            drawnow;
        end
        %%%%%%%%%%%%%
    end
    

    cost_red = 100*(cost_prev-cost)/cost_prev; % cost reduction in percent
    % print cost if demanded
    if verbose
        fprintf('Current cost: %.5e, stepsize: %.3e, iteration: %3.0f, cost reduction: %3.4f%%, runtime: %.1f\n', cost, s, iter-1, cost_red,toc);
    end
end


mrxisetup_optim = mrxisetup.copySetup;
cost_optim = cost;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print final cost if demanded
if verbose
    fprintf('Final cost: %.5e, runtime: %.1f\n', cost_optim, toc);
end

%%%%%%%%%%%%% debugging
if debug
    figure(121);
    scatter3(mrxisetup.coilData.centers(1,1),...
            mrxisetup.coilData.centers(2,1),...
            mrxisetup.coilData.centers(3,1),'filled','g');
    quiver3(mrxisetup.coilData.centers(1,1),...
        mrxisetup.coilData.centers(2,1),...
        mrxisetup.coilData.centers(3,1),...
        mrxisetup.coilData.orientations(1,1)*quiverScale,...
        mrxisetup.coilData.orientations(2,1)*quiverScale,...
        mrxisetup.coilData.orientations(3,1)*quiverScale,'g');

    figure(122);
    scatter3(mrxisetup.coilData.centers(1,1),...
            mrxisetup.coilData.centers(2,1),...
            mrxisetup.coilData.centers(3,1),'filled','g');
    quiver3(mrxisetup.coilData.centers(1,1),...
        mrxisetup.coilData.centers(2,1),...
        mrxisetup.coilData.centers(3,1),...
        mrxisetup.coilData.orientations(1,1)*quiverScale,...
        mrxisetup.coilData.orientations(2,1)*quiverScale,...
        mrxisetup.coilData.orientations(3,1)*quiverScale,'g');
    
    figure(123);
    scatter3(heightsCurr, anglesCurr, log10(cost), 'filled', 'g');
    drawnow;
end
%%%%%%%%%%%%%

end

function [L, L_notNormed, froL_single] = getLeadfieldMatrixWithNewCoilPositions(mrxisetup_in, VS, centers, orientations, normL_flag)
    % the function name says it all
    mrxisetup = mrxisetup_in.copySetup;
    
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
    
    mrxisetup.setCoils(centers, orientations, ...
        mrxisetup.coilData.coilpattern, ...
        mrxisetup.coilData.coilpatternassignment);
    mrxisetup.setCoilRadii(mrxisetup_in.getCoilRadii);
    
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


function ccenter_grad = updatePositionGradient(...
    mrxisetup, mrxisetup_init, constraintData, angles_init, angles_curr,...
    VS, L, pL, froL, fropL, froL_single, L_notNormed)
    % returns the derivative of the Frobenius condition of L w.r.t. the
    % coil positions on a cylindrical constraint 
    % OUTPUT: nCoils vector with gradient of cond(L,'fro') w.r.t. the coil
    % positions

    nCoils = mrxisetup.getNumberOfCoils;
    nSensors = mrxisetup.getNumberOfSensors;
    nVoxels = mrxisetup.getNumberOfVoxels;
    normL_flag = ~isempty(froL_single);
    
    p = permute(mrxisetup.ROIData.voxels,[1 3 2]);
    
    ccenter_grad = zeros(2,nCoils);
    
    pLpLT = pL*pL';
    npL_div_nL = fropL/froL;
    nL_div_npL = 1/npL_div_nL;
    
    % get derivative of rotation matrix w.r.t. angle for further
    % computation
    dRotMatrices = getDerivativeOfRotationMatricesFromInitPosition(...
        angles_curr, angles_init, constraintData);

    for ci = 1:nCoils
        % get current plane vectors for the right constraint
        coIdx = constraintData.constraintAlloc(ci);
        coDataIdx = find(constraintData.usedConstraints==coIdx);
        cpIdx = mrxisetup.coilData.coilpatternassignment(ci);
        nSegments = size(mrxisetup.coilData.coilpattern{cpIdx},2);
        % get segment vectors s of current coil
        s1 = bsxfun(@minus, mrxisetup.coilData.segments.(['c' num2str(ci)]), p);
        s2 = circshift(s1, -1, 2);
        % delete last segment points because they would connect the last
        % point to the first
        
        % compute derivatives of coil segments wrt. angle shifts (derivatives
        % wrt. height are simply the normalized constrained axes for each
        % segment)
        % derivative of coil segments w.r.t. height
        sder_height = normcols(constraintData.axesConst(:,coDataIdx));
        % derivative of coil segments w.r.t. angle
        c_init = mrxisetup_init.coilData.rotMatrices(:,:,ci) * ...
            mrxisetup_init.coilData.coilpattern{cpIdx}; % rotated coil pattern
        sder_angle1 = repmat(constraintData.rotMatrixConst(:,:,coDataIdx) * ...
            mrxisetup.coilData.constraints{coIdx}.radius * ...
            [-sin(angles_curr(ci)) cos(angles_curr(ci)) 0]',1,nSegments) + ...
            dRotMatrices(:,:,ci) * c_init;
        %%% multiply angle derivative with 1/constraintRadius^2 (I cannot
        %%% figure out why 1/constraintRadius^2 and not
        %%% 1/constraintRadius...)
        sder_angle1 = repmat(sder_angle1,1,1,nVoxels);
        sder_angle2 = circshift(sder_angle1, -1, 2);

        % compute norms of s
        ns1 = sqrt(sum(s1.*s1, 1));
        ns2 = circshift(ns1, -1, 2);
        % compute sum of norms
        ns1_p_ns2 = ns1 + ns2;
        % compute product of norms ||s1||*||s2||
        ns1ns2 = ns1.*ns2;
        % compute s/ns and then s1/ns1 + s2/ns2
        s1_d_ns1 = bsxfun(@rdivide, s1, ns1);
        % compute cross products s1 x s2
        CP_s1_s2 = cross(s1, s2, 1);
        % compute scalar products <s1, s2>
        SP_s1_s2 = sum(s1 .* s2, 1);

        %%% compute terms of the derivative of H step by step
        % numerator
        Num =   bsxfun(@times, ns1_p_ns2, CP_s1_s2);
        % denominator
        Den =   ns1ns2 .* (ns1ns2 + SP_s1_s2);

        for m = 1:2
            switch m
                case 1 % for height derivative
                    % get right size of sder_height
                    sder_height_matrix = repmat(sder_height, 1, nSegments, nVoxels);
                    % compute scalar products <s/||s||, ds>
                    SP_s1dns1_ds1 = sum(s1_d_ns1 .* sder_height_matrix, 1);
                    SP_s2dns2_ds2 = circshift(SP_s1dns1_ds1, -1, 2);
                    % derivative of numerator
                    dNum =  bsxfun(@times, SP_s1dns1_ds1 + SP_s2dns2_ds2, CP_s1_s2) + ...
                            bsxfun(@times, ns1_p_ns2, cross(sder_height_matrix, s2 - s1, 1));
                    % derivative of denominator
                    dDen =  (2*ns1ns2 + SP_s1_s2).*(SP_s1dns1_ds1.*ns2 + SP_s2dns2_ds2.*ns1) + ...
                            ns1ns2 .* sum(sder_height_matrix .* (s1 + s2), 1); 
                case 2 % for angle derivative
                    % compute scalar products <s/||s||, ds>
                    SP_s1dns1_ds1 = sum(s1_d_ns1 .* sder_angle1, 1);
                    SP_s2dns2_ds2 = sum(circshift(s1_d_ns1,-1,2) .* sder_angle2, 1);
                    % derivative of numerator
                    dNum =  bsxfun(@times, SP_s1dns1_ds1 + SP_s2dns2_ds2, CP_s1_s2) + ...
                            bsxfun(@times, ns1_p_ns2, cross(sder_angle1, s2, 1) + cross(s1, sder_angle2, 1));
                    % derivative of denominator
                    dDen =  (2*ns1ns2 + SP_s1_s2).*(SP_s1dns1_ds1.*ns2 + SP_s2dns2_ds2.*ns1) + ...
                            ns1ns2 .* (sum(sder_angle1 .* s2, 1) + sum(sder_angle2 .* s1, 1)); 
            end
            

            %%% compute gradient of H
            H_grad =...
                1/(4*pi)*sum(bsxfun(@rdivide,...
                (bsxfun(@times, dNum(:,1:end-1,:), Den(:,1:end-1,:)) - bsxfun(@times, Num(:,1:end-1,:), dDen(:,1:end-1,:))),...
                Den(:,1:end-1,:).^2),2); 

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
            ccenter_grad(m,ci) = sum(sum(L(ii,:).*L_grad))*npL_div_nL - nL_div_npL*sum(sum(pLpLT.*(pL(:,ii)*L_grad)));
            if m == 2
                ccenter_grad(m,ci) = ccenter_grad(m,ci)./mrxisetup.coilData.constraints{coIdx}.radius^2;
            end
        end
    end
end


function [anglesOnConst, heightOnCons] = getCylindricalCoords(mrxisetup, constraintData)
% returns the cylindrical coordinates on constraints

% initialize variables
nCoils = length(constraintData.constraintAlloc);
anglesOnConst = zeros(nCoils,1);
heightOnCons = zeros(nCoils,1);

% for each coil return angle and height on constraint
for c = 1:nCoils
    % get index of currently used constraint
    ci = constraintData.constraintAlloc(c);
    coIdx = find(constraintData.usedConstraints == ci);
    % get current coil center
    ccenter = mrxisetup.coilData.centers(:,c);
    % get constraint bottom coordinate and axis vectors
    a0 = mrxisetup.coilData.constraints{ci}.bottomCoord;
    coAxisNormed = normcols(constraintData.axesConst(:,coIdx));
    coRotNormed = constraintData.angleRefAxesConst(:,coIdx);
    
    % compute height on cylinder axis
    heightOnCons(c) = coAxisNormed' * (ccenter - a0);
    
    % get projected coil center on constraint axis
    projCenter = a0 + heightOnCons(c)*coAxisNormed;
    % get normalized vector from projected point to actual coil center
    projToCenter = normcols(ccenter - projCenter);
    % get positive angle direction
    coDirNormed = constraintData.angleDirAxesConst(:,coIdx);
    % compute angle
    angleOneToPi = acos(dot(coRotNormed, projToCenter));
    % add Pi if necessary
    if dot(coDirNormed, projToCenter) < 0
        anglesOnConst(c) = 2*pi - angleOneToPi;
    else
        anglesOnConst(c) = angleOneToPi;
    end
end
end

function [centers, orientations] = getEuclideanCoordinates(...
    mrxisetup_init, anglesOnConst_init, anglesOnConst, heightOnCons, constraintData)
% returns the euclidean coordinates

% initialize variables
nCoils = length(constraintData.constraintAlloc);
centers = zeros(3,nCoils);
orientations = zeros(3,nCoils);

% compute rotation matrices for orientations
rotMatricesFromInitToCurr = getRotationMatricesFromInitPosition(...
    anglesOnConst, anglesOnConst_init, constraintData);


% for each coil return angle and height on constraint
for c = 1:nCoils
    % get index of currently used constraint
    ci = constraintData.constraintAlloc(c);
    coIdx = constraintData.usedConstraints == constraintData.constraintAlloc(c);
    % get constraint rotation matrix
    rotMatrixConst = constraintData.rotMatrixConst(:,:,coIdx);
    % get bottom coordinate of cylindrical constraint
    a0 = mrxisetup_init.coilData.constraints{ci}.bottomCoord;
    % get cylinder radius
    r = mrxisetup_init.coilData.constraints{ci}.radius;
    
    % compute new center
    centers(:,c) = a0 + ...
        rotMatrixConst*[r*cos(anglesOnConst(c)); r*sin(anglesOnConst(c)); heightOnCons(c)];
    
    % compute new orientations
    orientations(:,c) = rotMatricesFromInitToCurr(:,:,c) * ...
        mrxisetup_init.coilData.orientations(:,c);
end
end


function rotMatrices = getRotationMatricesFromInitPosition(angles, anglesInit, constraintData)
% returns the rotation matrices to rotate from initial positions of coils
% to position defined by angles; also returns the deviation of the rotation
% matrices with respect to the angles if demanded

% initialize variables
nCoils = length(constraintData.constraintAlloc);
rotMatrices = zeros(3,3,nCoils);

% for each coil return angle and height on constraint
for c = 1:nCoils
    % get index of currently used constraint
    coIdx = constraintData.usedConstraints == constraintData.constraintAlloc(c);
    % get current normalized constraint axis
    cAx = normcols(constraintData.axesConst(:,coIdx));
    % get angle deviation from initial position
    da = angles(c) - anglesInit(c);
    % calculate rotation matrix
    cosda = cos(da);
    sinda = sin(da);
    omcosda = 1 - cosda;
    rotMatrices(:,:,c) = ...
        [cAx(1)^2*omcosda+cosda  cAx(1)*cAx(2)*omcosda-cAx(3)*sinda  cAx(1)*cAx(3)*omcosda+cAx(2)*sinda; ...
         cAx(1)*cAx(2)*omcosda+cAx(3)*sinda  cAx(2)^2*omcosda+cosda  cAx(2)*cAx(3)*omcosda-cAx(1)*sinda; ...
         cAx(1)*cAx(3)*omcosda-cAx(2)*sinda  cAx(2)*cAx(3)*omcosda+cAx(1)*sinda  cAx(3)^2*omcosda+cosda];
end


end


function dRotMatrices = getDerivativeOfRotationMatricesFromInitPosition(angles, anglesInit, constraintData)
% returns the rotation matrices to rotate from initial positions of coils
% to position defined by angles; also returns the deviation of the rotation
% matrices with respect to the angles if demanded

% initialize variables
nCoils = length(constraintData.constraintAlloc);
dRotMatrices = zeros(3,3,nCoils);

% for each coil return angle and height on constraint
for c = 1:nCoils
    % get index of currently used constraint
    coIdx = constraintData.usedConstraints == constraintData.constraintAlloc(c);
    % get current normalized constraint axis
    cAx = normcols(constraintData.axesConst(:,coIdx));
    % get angle deviation from initial position
    da = angles(c) - anglesInit(c);
    % calculate rotation matrix
    cosda = cos(da);
    sinda = sin(da);
    dRotMatrices(:,:,c) = ...
        [cAx(1)^2*sinda-sinda  cAx(1)*cAx(2)*sinda-cAx(3)*cosda  cAx(1)*cAx(3)*sinda+cAx(2)*cosda; ...
         cAx(1)*cAx(2)*sinda+cAx(3)*cosda  cAx(2)^2*sinda-sinda  cAx(2)*cAx(3)*sinda-cAx(1)*cosda; ...
         cAx(1)*cAx(3)*sinda-cAx(2)*cosda  cAx(2)*cAx(3)*sinda+cAx(1)*cosda  cAx(3)^2*sinda-sinda];
end


end




