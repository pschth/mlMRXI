function [SETUP, optimizationPath] = optimizeCylindricalCoilPositionsFroCond( SETUP, InitialPositions, s_min, s_max, s_interSteps, maxIter, verbose, debug )
% optimizes the coil positions using the Frobenius condition of the
% Leadfield matrix as objective function to be minimized. 
% 
% INPUT:
% SETUP - MRXI setup configuration object
% 
% InitialPositions - 2 x nCoils matrix; initial coil positions specified by
% angles (radians) in the first matrix row and the distances along the vector from
% the bottom to the top of the cylinder
% 
% s_min (optional) ... minimum stepsize, default: 1e-3
% 
% s_max (optional) ... maximum stepsize, default: 1
% 
% s_interSteps (optional) ... number of logarithmically equidistant steps
% between s_min and s_max, default: 25
% 
% maxIter (optional) ... maximum number of iterations, default: Inf 
% 
% verbose (optional) - boolean; prints optimization step information;
% default: true
% 
% debug (optional) - boolean; plots the setup and the optimization path;
% default: false

% exceptions & default values
if nargin < 3
    s_min = 1e-3;
elseif isempty(s_min)
    s_min = 1e-3;
end

if nargin < 4
    s_max = 1e-0;
elseif isempty(s_max)
    s_max = 1e-0;
end

if nargin < 5
    s_interSteps = 25;
elseif isempty(s_interSteps)
    s_interSteps = 25;
end

if nargin < 6
    maxIter = Inf;
elseif isempty(maxIter)
    maxIter = Inf;
end

if nargin < 7
    verbose = true;
elseif isempty(verbose)
    verbose = true;
end

if nargin < 8
    debug = false;
elseif isempty(debug)
    debug = false;
end


% % check if there is a (only one) cylindrical constraint specified on the
% SETUP
if ~isfield(SETUP.coilData, 'constraints')
    error('No cylindrical constraint on SETUP specified');
else
    cylindricalConstraintIdx = false(length(SETUP.coilData.constraints),1);
    for i = 1:length(SETUP.coilData.constraints)
        if strcmp(SETUP.coilData.constraints{i}.type, 'cylindrical')
            cylindricalConstraintIdx(i) = true;
        end
    end
    if sum(cylindricalConstraintIdx)~= 1
        error('No or more than one cylindrical constraint found. Only one allowed.');
    else
        % get the valid constraint
        constraint = SETUP.coilData.constraints{cylindricalConstraintIdx};
    end
end

% check right size of InitialPositions
if ~any(size(InitialPositions) == 2)
    error('InitialPositions does not have 2 rows (or columns).');
elseif size(InitialPositions,1) ~= 2 && size(InitialPositions,2) == 2
    % transpose if necessary
    InitialPositions = InitialPositions';
end
positions = InitialPositions;

% check verbose
if nargin < 5
    verbose = true;
elseif isempty(verbose)
    verbose = true;
end

% check debug
if nargin < 6
    debug = false;
elseif isempty(debug)
    debug = false;
end

%%% code
tic % start timer
nCoils = size(positions,2);
nSensors = SETUP.getNumberOfSensors;

% store each optimization step if required
if nargout == 2
    storePath = true;
    optimizationPath = zeros(2,nCoils,200);
    optimizationPath(:,:,1) = positions;
else
    storePath = false;
end

% get the cylinder axis vector and translation offset
cAxis = constraint.topCoord - constraint.bottomCoord;

% check if InitialPositions are off the constraint
axisLength = norm(cAxis);
for c = 1:nCoils
    if positions(2,c) > axisLength
        error(['InitialPosition no. ' num2str(c) ' off the constraint. Maximum bottom-top distance: ' num2str(axisLength)]);
    end    
end

% set the InitialPositions in the SETUP object
    % function to set points on cylindrical constraint
[centers, orientations, rotMatrix] = getEuclideanCoords(positions, constraint);
    % set coils
if SETUP.coilData.isDipole
    SETUP.setCoils(centers, orientations);
else
    coilpattern = SETUP.coilData.coilpattern;
    coilpatternassignment = SETUP.coilData.coilpatternassignment;
    SETUP.setCoils(centers, orientations, coilpattern, coilpatternassignment);
end

%%%%%%%%%%%%%%% for debugging
if debug
    fig1 = figure(501); clf;
    SETUP.visualize(true);
    hold on;
    visualizeCostFunctionOfOneCoil( SETUP, [], 25 );
    visualizeCoilGrid(SETUP, 'Color', 'm');
    pause(0.1)
    fig2 = figure(502); clf;
    hold on;
    visualizeCostFunctionOfOneCoil( SETUP, 'flat', 25 );
    set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi])
    set(gca,'XTickLabel',{'0' '\pi/2' '\pi' '3\pi/2' '2\pi'})
    title('Frobenius condition number ||L||_F||L^{-1}||_F');
    scatter(mod(positions(1,:), 2*pi), positions(2,:), 'filled', 'm');
    pause(0.1);
    figure(503); clf;
    visualizeCostFunctionOfOneCoil( SETUP, 'flat', 25, false );
    title('Spectral condition number ||L||_2||L^{-1}||_2');
    set(gca,'XTick',[0 pi/2 pi 3*pi/2 2*pi])
    set(gca,'XTickLabel',{'0' '\pi/2' '\pi' '3\pi/2' '2\pi'})
    pause(0.1);
end
%%%%%%%%%%%%%%%

% get voxel-sensor relationship matrix (static)
S = getVoxelSensorRelation(SETUP,false);

% set sequential activation pattern
I = speye(nCoils);

% compute initial leadfield matrix
L = createSystemMatrix( SETUP, I, [], false );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimization routine

% initialize
iter = 1; % counter
cost = algoCost( L ); % cost of previous step
if verbose
    fprintf('Initial cost: %.4e\n', cost);
end

grad = zeros(2,nCoils); % actual gradient

s_search = logspace(log10(s_min), log10(s_max), s_interSteps);
cost_temp = zeros(s_interSteps, 1);

while true
    iter = iter + 1; % increment counter
    if iter > maxIter
        if verbose
            fprintf('Hit maximum number of iterations %.0f. Terminating optimization.\n', maxIter);
        end
        break
    end
    
%     % get voxel coil relationship matrix
%     C = getVoxelCoilRelation(SETUP, false);
    % get partial derivative matrices w.r.t. angle and distance of voxel
    % coil relationship matrix
    C_dot = getGradientOfVoxelCoilRelation(SETUP, positions, constraint, rotMatrix);
    L_dot = getGradientOfLeadfieldMatrix(S, C_dot);
    
    % get pseudoinverse of current L
    pL = pinv(L);
    
    % get Frobenius norms of L and pL
    froL = norm(L,'fro');
    fropL = norm(pL,'fro');

    % calculation of the actual gradient
    
    for c = 1:nCoils
        blockIdx = (c-1)*nSensors+1:c*nSensors;
        for a = 1:2
            leftTerm = sum(sum(L(blockIdx,:).*L_dot(:,:,c,a)))*fropL/froL;
            rightTerm = trace(pL'*(pL(:,blockIdx)*L_dot(:,:,c,a))*pL)*froL/fropL;
            grad(a,c) = leftTerm - rightTerm;
        end
    end
    
    % update positions
        % multiply factor for cylindricial coordinate derivative
        %%%%%%%%%%% NOTE:
        % factor 1/Radius of the angle derivative is mathematically correct
        % but a re-conversion into radians (factor 2*Pi/Radius) works much
        % better and converges much faster! 
        %%%%%%%%%%%
%     grad(1,:) = grad(1,:)*1./constraint.radius; % correct version
    grad(1,:) = grad(1,:)*2*pi./constraint.radius; % altered version
    grad = grad./norm(grad);
    
    % determine stepsize with largest cost reduction within s_min and s_max
    % using s_interSteps with logarithmically equal 
    for si = 1:s_interSteps
        positions_temp = positions - s_search(si)*grad;
        % set new positions in SETUP object
        [centers, orientations] = getEuclideanCoords(positions_temp, constraint, rotMatrix);
        if SETUP.coilData.isDipole
            SETUP.setCoils(centers, orientations);
        else
            coilpattern = SETUP.coilData.coilpattern;
            coilpatternassignment = SETUP.coilData.coilpatternassignment;
            SETUP.setCoils(centers, orientations, coilpattern, coilpatternassignment);
        end
        % update cost
        cost_temp(si) = algoCost( createSystemMatrix( SETUP, I, [], false ) );
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
        positions = positions - s*grad;
        % set new coil positions
        [centers, orientations, rotMatrix] = getEuclideanCoords(positions, constraint);
        if SETUP.coilData.isDipole
            SETUP.setCoils(centers, orientations);
        else
            coilpattern = SETUP.coilData.coilpattern;
            coilpatternassignment = SETUP.coilData.coilpatternassignment;
            SETUP.setCoils(centers, orientations, coilpattern, coilpatternassignment);
        end
        % calculate new leadfield matrix
        L = createSystemMatrix( SETUP, I, [], false );
        
        if storePath
            optimizationPath(:,:,iter) = positions;
        end
        cost_prev = cost;
        cost = cost_temp(si-1);
    end
    
    % print cost if demanded
    if verbose
        fprintf('Current cost: %.4e, stepsize: %.3e, iterations: %.0f, cost reduction: %3.4f%%, runtime: %.1f\n', cost, s, iter-1, 100*(cost_prev-cost)/cost_prev, toc);
    end
    
    %%%%%%%%%%%%%%% for debugging
    if debug
        figure(fig1)
        hold on
        visualizeCoilGrid(SETUP, 'Color', 'r');
        drawnow
        figure(fig2)
        scatter(mod(positions(1,:), 2*pi), positions(2,:), 'filled', 'r');
        drawnow
    end
    %%%%%%%%%%%%%%%
end

% print final cost if demanded
if verbose
    fprintf('Final cost: %.4e, iterations: %.0f, runtime: %.1f\n', cost, iter-1, toc);
end

%%%%%%%%%%%%%%% for debugging
if debug
    figure(fig1)
    hold on
    visualizeCoilGrid(SETUP, 'Color', 'k');
    pause(0.01)
    figure(fig2)
    scatter(mod(positions(1,:), 2*pi), positions(2,:), 'filled', 'k');
end
%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% other functions
function cost = algoCost( L, tol )
    if nargin == 4
        L = L./normest(L);
        pL = pinv( L ,tol);
    else
        pL = pinv(L);
    end
    cost = norm(L,'fro') * norm(pL,'fro');
end


function C_dot = getGradientOfVoxelCoilRelation(SETUP, positions, constraint, rotMatrix)
% returns the partial derivatives of the voxel coil relationship matrix
% with respect to the cylindrical constraints (i.e. the partial derivatives
% w.r.t. the angle and the distance along the cylinder axis)
% 
% INPUT:
% SETUP - MRXI setup configuration object (MUST already contain the right
% coil positions on the cylindrical constraint)
% 
% Positions - 2 x nCoils matrix; coil positions specified by
% angles (radians) in the first matrix row and the distances along the vector from
% the bottom to the top of the cylinder
% 
% rotMatrix - 3 x 3 matrix; rotation matrix from [0 0 1] to the actual axis
% orientation of the cylindrical constraint
% 
% OUTPUT:
% C_dot - 3 x nCoils x nVoxels x 2 matrix; the two 3 x nCoils x nVoxels
% gradient matrices of the voxel coil relationship w.r.t. (1) angle and (2)
% distance along the cylinder axis

%%%%%%%%%%% NOTE:
% rc = SETUP.coilData.centers; coil vectors
% nc = SETUP.coilData.orientations; coil orientation vectors
% rk = SETUP.ROIData.voxels; voxel vectors
%%%%%%%%%%%

% get some constants
nCoils = SETUP.getNumberOfCoils;
nVox = SETUP.getNumberOfVoxels;

% initialize C_dot
C_dot = zeros(3,nCoils,nVox,2);

% calculate partial derivatives of coil positions and orientations
dr0 = [ -constraint.radius*sin(positions(1,:));...
        constraint.radius*cos(positions(1,:));...
        zeros(1,nCoils)]; % PD of unrotated coil positions (along [0 0 1])
drca = - rotMatrix * dr0; % partial derivative of (rk-rc) w.r.t. angle
dnca = drca./constraint.radius; % partial derivative of nc w.r.t. angle
drcz = - rotMatrix * repmat([0;0;1],1,nCoils); % partial derivative of (rk-rc) w.r.t. distance
%%%%%%%%%%% NOTE:
% dncz = 0; % partial derivative of nc w.r.t. distance is zero
%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% calculation of partial derivative
for c = 1:nCoils
    % get required coil vectors
    rc = SETUP.coilData.centers(:,c);
    nc = SETUP.coilData.orientations(:,c);
    darkc = drca(:,c);
    dankc = dnca(:,c);
    dzrkc = drcz(:,c);
    for k = 1:nVox
        rk = SETUP.ROIData.voxels(:,k);
        rkc = rk-rc;
        rkc_norm = norm(rkc);
        
        % PDE w.r.t. angle
        C_dot(:,c,k,1) = 1/(4*pi) * (...
            3*((darkc*nc'*rkc+rkc*dankc'*rkc+rkc*nc'*darkc)*rkc_norm^2 - ...
            5*(darkc'*rkc)*rkc*nc'*rkc)./rkc_norm^7 -...
            (dankc*rkc_norm^2 - 3*(darkc'*rkc)*nc)./rkc_norm^5 );
        % PDE w.r.t. distance
        C_dot(:,c,k,2) = 3/(4*pi) * (...
            ((dzrkc*nc'*rkc+rkc*nc'*dzrkc)*rkc_norm^2 - ...
            5*(dzrkc'*rkc)*rkc*nc'*rkc)./rkc_norm^7 +...
            ((dzrkc'*rkc)*nc)./rkc_norm^5 );
    end
end

end


function L_dot = getGradientOfLeadfieldMatrix(S, C_dot)
% returns the partial derivatives of the leadfield matrix with respect to
% the cylindrical constraints (i.e. the partial derivatives w.r.t. the
% angle and the distance along the cylinder axis) 
% 
% INPUT:
% S - 3 x nSensors x nVoxel matrix; voxel-sensor relationship matrix
% 
% C_dot - 3 x nCoils x nVoxel x 2 matrix; derived voxel-coil relationship
% matrices w.r.t. angle and distance along the cylinder axis
% 
% OUTPUT:
% L_dot: nSensors x nVoxel x nCoils x 2 matrix; returns the derived
% leadfield matrix w.r.t. all coil positions which depend on angle and
% distance along the cylinder axis, respectively (=nCoils x 2 derivations)

[~, nSensors, nVoxel] = size(S);
nCoils = size(C_dot,2);

L_dot = zeros(nSensors, nVoxel, nCoils, 2);

for c = 1:nCoils
    L_dot(:,:,c,:) = permute(sum(bsxfun(@times, S, C_dot(:,c,:,:)),1), [2 3 1 4]);
end

end

function [centers, orientations, rotMatrix] = getEuclideanCoords(positions, constraint, rotMatrix)
% calculates the coil positions (centers) and their orientations in
% euclidean coordinates from cylinder coordinates (positions(1,:) = angles,
% positions(2,:) = z-distance).

r0 =  [ constraint.radius*cos(positions(1,:));...
        constraint.radius*sin(positions(1,:));...
        positions(2,:)]; % cylinder coordinates to euclidean on cylinder along [0 0 1] axis

d0 = [r0(1:2,:); zeros(1, size(positions,2)) ]; % coil orientation vectors on cylinder along [0 0 1] axis

% rotate points on cylinder
if nargin < 3
    direction = normcols(constraint.topCoord - constraint.bottomCoord);
    [centers, rotMatrix] = rotatepoints(r0, direction);
else
    centers = rotMatrix*r0;
end

% translate to origin of constraint
centers = bsxfun(@plus, centers, constraint.bottomCoord);
% rotate orientations
orientations = -rotMatrix*d0;

end


function [angleGridpoints, zGridpoints, gx, gy] = visualizeCostFunctionOfOneCoil( SETUP, representation, nSteps, cost_or_cond )
% visualizes the cost function of a SETUP with a single dipole coil along a
% cylindrical constraint
% 
% INPUT
% SETUP - MRXI setup configuration object
% 
% representation (optional) - string; either 'flat' or 'cylindrical',
% defines the representation of the cost distribution. 'flat' makes a 2D
% plot with an angle- and a z-axis. 'cylindrical' is a 3D cylinder;
% default: 'cylindrical'
% 
% nSteps(optional) - scalar or 2 x 1 vector; defines the number of steps in
% (1) angle direction and (2) axis direction
% 
% cost_or_cond (optional) - boolean; true: Frobenius condition, false: spectral
% condition, Default: true


% check if there is a (only one) cylindrical constraint specified on the
% SETUP
if ~isfield(SETUP.coilData, 'constraints')
    error('No cylindrical constraint on SETUP specified');
else
    cylindricalConstraintIdx = false(length(SETUP.coilData.constraints),1);
    for i = 1:length(SETUP.coilData.constraints)
        if strcmp(SETUP.coilData.constraints{i}.type, 'cylindrical')
            cylindricalConstraintIdx(i) = true;
        end
    end
    if sum(cylindricalConstraintIdx)~= 1
        error('No or more than one cylindrical constraint found. Only one allowed.');
    else
        % get the valid constraint
        constraint = SETUP.coilData.constraints{cylindricalConstraintIdx};
    end
end

if nargin < 4
    cost_or_cond = true;
end

if nargin < 2
    representation = 'cylindrical';
elseif isempty(representation)
    representation = 'cylindrical';
elseif ~(strcmp(representation,'cylindrical') || strcmp(representation,'flat'))
    error('''representation'' has to be either ''flat'' or ''cylindrical''.');
end

% copy SETUP to not change the settings
mrxisetup = SETUP.copySetup();
nCoils = size(mrxisetup.coilData.centers,2);

% get the cylinder axis vector and translation offset
cAxis = constraint.topCoord - constraint.bottomCoord;
cOffset = constraint.bottomCoord;
% get the normalized cylinder axis to rotate to
direction = normcols(cAxis);
% get z length
zl = norm(cAxis);


%% SETUP CONFIGURATION
% set visualization grid properties
if nargin < 2
    angleSteps = 20;
    zSteps = 20;
elseif isempty(nSteps)
    angleSteps = 20;
    zSteps = 20;
elseif numel(nSteps)==1
    angleSteps = nSteps;
    zSteps = nSteps;
elseif numel(nSteps)==2
    angleSteps = nSteps(1);
    zSteps = nSteps(2);
else
    error('Wrong input size of nSteps. Either scalar or 2 x 1 vector.');
end
angleGridpoints = linspace(0,2*pi,angleSteps);
zGridpoints = linspace(0, zl, zSteps);
[angleGrid, zGrid] = meshgrid(angleGridpoints, zGridpoints);

if nCoils == 1
    % set initial positions
    positions = [angleGrid(:)'; zGrid(:)'];
    nGridpoints = size(positions,2);

    % initialize cost calculation loop
    centers = zeros(3, nGridpoints);
    cost = zeros(1, nGridpoints);
    % cost calculation loop
    for p = 1:size(positions,2)
        % set points on circle oriented along the z-axis
        cylPositions = @(positions) [   constraint.radius*cos(positions(1,:));...
                                        constraint.radius*sin(positions(1,:));...
                                        0];
        centers(:,p) = cylPositions(positions(:,p));
        % the coils are oriented towards the center of the cylinder,
        % perpendicular to the cylinder surface
        cylOrientations = @(centers) -normcols(centers);
        orientations = cylOrientations(centers(:,p));
        % translate points along z-axis
        centers(:,p) = centers(:,p) + [0;0;1]*positions(2,p);
        % rotate the points and orientations and store the rotation matrix for
        % later use
        if p == 1
            [centers(:,p), rotMatrix] = rotatepoints(centers(:,p), direction);
            orientations = rotatepoints(orientations, direction);
        else
            centers(:,p) = rotMatrix*centers(:,p);
            orientations = rotMatrix*orientations;
        end
        % translate to origin of constraint
        centers(:,p) = bsxfun(@plus, centers(:,p), cOffset);
        if nCoils == 1
            % set coils
            mrxisetup.setCoils(centers(:,p), orientations);
            % calculate cost of new leadfield matrix
            if cost_or_cond
                cost(p) = algoCost( createSystemMatrix( mrxisetup, 1, [], false ) );
            else
                cost(p) = cond( createSystemMatrix( mrxisetup, 1, [], false ) );
            end
        end
    end
end

if nargout > 2 && nCoils == 1
    [gx,gy] = gradient(reshape(cost, size(angleGrid)));
else
    gx = [];
    gy = [];
end
    
    
if strcmp(representation, 'cylindrical')
    %% visualize cylindrical cost distribution of cylindrical coil positioning
    if nCoils == 1
        X = reshape(centers(1,:), size(angleGrid));
        Y = reshape(centers(2,:), size(angleGrid));
        Z = reshape(centers(3,:), size(angleGrid));
        cost = reshape(cost, size(angleGrid));
        surf(X,Y,Z,reshape(cost, size(X)));
        shading interp;
        colorbar;
    else
        mrxisetup.visualizeCoilConstraints(50);
    end
    % set(p1,'FaceAlpha',0.7);
    daspect([1 1 1]);
    view([-45 25]);
    xlabel('x'); ylabel('y'); zlabel('z');
else
    if nCoils == 1
        pcolor(angleGrid, zGrid, reshape(cost, size(zGrid)));
        shading interp;
        colorbar;
    else
        axis;
    end
    xlabel('angle (rad)'); ylabel('z-distance');
    xlim([0 2*pi]);
    ylim([0 zl]);
end
drawnow;
end

