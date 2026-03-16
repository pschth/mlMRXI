function [SETUP, optimizationPath] = optimizeCylindricalCoilPositionsXYZ( SETUP, InitialPositions, s_init, s_thresh, c_armijo, verbose )
% optimizes the coil positions using the Frobenius norm of the Gram matrix
% of the Leadfield matrix as cost function to be minimized.
% 
% INPUT:
% SETUP - MRXI setup configuration object
% 
% InitialPositions - 3 x nCoils matrix; initial coil positions
% 
% s_init ... initial gradient-descent stepsize (recommendation: 1e-1)
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
if ~any(size(InitialPositions) == 3)
    error('InitialPositions does not have 2 rows (or columns).');
elseif size(InitialPositions,1) ~= 3 && size(InitialPositions,2) == 3
    % transpose if necessary
    InitialPositions = InitialPositions';
end
positions = InitialPositions;

% check verbose
if nargin < 6
    verbose = true;
elseif isempty(verbose)
    verbose = true;
end

%%% code
tic % start timer
nCoils = size(positions,2);
nVox = SETUP.getNumberOfVoxels;

% store each optimization step if required
if nargout == 2
    storePath = true;
    optimizationPath = zeros(3,nCoils,200);
    optimizationPath(:,:,1) = positions;
else
    storePath = false;
end

% get the cylinder axis vector and translation offset
cAxis = constraint.topCoord - constraint.bottomCoord;
ncAxis = cAxis./norm(cAxis);

% project all InitialPositions on the cylinder constraint
[centers, orientations] = projectOnCylinderSurface(positions, constraint, ncAxis);
% set coils
SETUP.setCoils(centers, orientations);

%%%%%%%%%%%%%%% for debugging
figure(1); clf;
SETUP.visualize(true);
hold on;
visualizeCostFunctionOfOneCoil( SETUP );
SETUP.setCoils(centers, orientations);
pause(0.1)
%%%%%%%%%%%%%%%

% get voxel-sensor relationship matrix (static)
S = getVoxelSensorRelation(SETUP,false);

% set sequential activation pattern
I = speye(nCoils);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimization routine

% initialize
s = s_init; % stepsize
cost_prev = froGramCost( createSystemMatrix( SETUP, I, [], false ) ); % cost of previous step
cost = cost_prev; % cost of current step
if verbose
    fprintf('Initial cost: %.6f\n', cost_prev);
end

ii = 1; % counter
while s > s_thresh
    ii = ii+1; % increment counter
    
    % get voxel coil relationship matrix
    C = getVoxelCoilRelation(SETUP, false);
    % get partial derivative matrices w.r.t. angle and distance of voxel
    % coil relationship matrix
    C_dot = getXYZGradientOfVoxelCoilRelation(SETUP);
    
    % precalculate all trace values required for gradient calculation
    % depending on i = 1:nVox
    Trii = zeros(nVox,1); % vector with all traces of Ci'*Si*Si'*Ci
    dTrii = zeros(3,nCoils,nVox); % vector with all traces of Ci_dot'*Si*Si'*Ci
    for i = 1:nVox% extract only angle gradient of the coil c
        Si = S(:,:,i)*S(:,:,i)';
        Trii(i) = trace(C(:,:,i)'*Si*C(:,:,i));
        for c = 1:nCoils
            % extract only gradient of the coil c
            for d = 1:3 % d ... dimension (angle & z-distance)
                dTrii(d,:,i) = diag(C_dot(:,:,i,d)'*Si*C(:,:,i));
            end
        end
    end
    
    % precalculate all trace values required for gradient calculation
    % depending on i = 1:nVox-1 and j=i:nVox
    Trij = zeros(nVox-1, nVox-1); % vector with all traces of Ci'*Si*Sj'*Cj
    Trdij = zeros(3, nCoils, nVox-1, nVox-1); % vector with all traces of Ci_dot'*Si*Sj'*Cj
    Tridj = Trdij; % vector with all traces of Ci'*Si*Sj'*Cj_dot
    for i = 1:nVox-1
        Ci = C(:,:,i);
        Si = S(:,:,i);
        for j = i+1:nVox
            Sij = Si*S(:,:,j)';
            CiSij = Ci'*Sij;
            SijCj = Sij*C(:,:,j);
            Trij(i,j-i) = trace(CiSij*C(:,:,j));
            % extract only angle gradient of the coil c
            for d = 1:3 % d ... dimension (angle & z-distance)
                Trdij(d,:,i,j-i) = diag(C_dot(:,:,i,d)'*SijCj);
                Tridj(d,:,i,j-i) = diag(CiSij*C_dot(:,:,j,d));
            end
        end
    end

    % calculation of the actual gradient
    grad = zeros(3,nCoils);
    for c = 1:nCoils
        ij = 0;
        for i = 1:nVox-1
            for j = i+1:nVox
                ij = ij+1;
                for d = 1:3 % d ... dimension (angle & z-distance)
                    % gradient calculation
                    grad(d,c) = grad(d,c) + ...
                        (Trij(i,j-i)*(Trdij(d,c,i,j-i) + Tridj(d,c,i,j-i)) - ...
                        (dTrii(d,c,i)/Trii(i) + dTrii(d,c,j)/Trii(j))*Trij(i,j-i)^2 )/...
                        (Trii(i)*Trii(j));
                end
            end
        end
    end
    
    % update positions
    grad = grad./norm(grad);
    positions_temp = positions - s*grad;
    
    % set new positions in SETUP object
    [centers, orientations] = projectOnCylinderSurface(positions_temp, constraint, ncAxis);
    SETUP.setCoils(centers,orientations);
    
    % update cost
    cost_prev = cost;
    cost = froGramCost( createSystemMatrix( SETUP, I, [], false ) );
    
    % print cost if demanded
    if verbose
        fprintf('Current cost: %.6f, armijo cost: %.6f, stepsize: %.3e, runtime: %.1f\n', cost, c_armijo, s, toc);
    end
    
    % stepsize control
    while cost_prev < cost
        s = s/2;
        if s < s_thresh
            break
        end
        % update positions
        positions_temp = positions - s*grad;
        [centers, orientations] = projectOnCylinderSurface(positions_temp, constraint, ncAxis); 
        % set new positions in SETUP object
        SETUP.setCoils(centers,orientations);
        
        % calculate cost of new leadfield matrix
        cost = froGramCost( createSystemMatrix( SETUP, I, [], false ) );
        % print cost if demanded
        if verbose
            fprintf('\tCurrent cost: %.6f, armijo cost: %.6f, stepsize: %.3e, runtime: %.1f\n', cost, c_armijo, s, toc);
        end
    end
    % update position and orientation if cost decreased
    positions = centers;
    
    % store optimization path if demanded
    if storePath
        optimizationPath(:,:,ii) = positions;
    end
    
    %%%%%%%%%%%%%%% for debugging
    figure(1);
    hold on
    visualizeCoilGrid(SETUP, 'Color', 'r');
    pause(0.01)
    %%%%%%%%%%%%%%%
end

end


function C_dot = getXYZGradientOfVoxelCoilRelation(SETUP)
% returns the partial derivatives of the voxel coil relationship matrix
% with respect to the cylindrical constraints (i.e. the partial derivatives
% w.r.t. the angle and the distance along the cylinder axis)
% 
% INPUT:
% SETUP - MRXI setup configuration object (MUST already contain the right
% coil positions on the cylindrical constraint)
% 
% OUTPUT:
% C_dot - 3 x nCoils x nVoxels x 3 matrix; the two 3 x nCoils x nVoxels
% gradient matrices of the voxel coil relationship w.r.t. (1) x, (2)
% y and (3) z

%%%%%%%%%%% NOTE:
% rc = SETUP.coilData.centers; coil vectors
% nc = SETUP.coilData.orientations; coil orientation vectors
% rk = SETUP.ROIData.voxels; voxel vectors
%%%%%%%%%%%

% get some constants
nCoils = SETUP.getNumberOfCoils;
nVox = SETUP.getNumberOfVoxels;

% initialize C_dot
C_dot = zeros(3,nCoils,nVox,3);

% calculate partial derivatives of coil positions and orientations
rx = -[1;0;0];
ry = -[0;1;0];
rz = -[0;0;1];
nx = rx;
ny = ry;
%%%%%%%%%%% NOTE:
% nz = 0;
%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% calculation of partial derivative
for c = 1:nCoils
    % get required coil vectors
    rc = SETUP.coilData.centers(:,c);
    nc = SETUP.coilData.orientations(:,c);
    for k = 1:nVox
        rk = SETUP.ROIData.voxels(:,k);
        rkc = rk-rc;
        rkc_norm = norm(rkc);
        
        % PDE w.r.t. x
        C_dot(:,c,k,1) = 1/(4*pi) * (...
            3*((rx*nc'*rkc+rkc*nx'*rkc+rkc*nc'*rx)*rkc_norm^2 - ...
            5*(rx'*rkc)*rkc*nc'*rkc)./rkc_norm^7 -...
            (nx*rkc_norm^2 - 3*(rx'*rkc)*nc)./rkc_norm^5 );
        % PDE w.r.t. y
        C_dot(:,c,k,2) = 1/(4*pi) * (...
            3*((ry*nc'*rkc+rkc*ny'*rkc+rkc*nc'*ry)*rkc_norm^2 - ...
            5*(ry'*rkc)*rkc*nc'*rkc)./rkc_norm^7 -...
            (ny*rkc_norm^2 - 3*(ry'*rkc)*nc)./rkc_norm^5 );
        % PDE w.r.t. z
        C_dot(:,c,k,3) = 3/(4*pi) * (...
            ((rz*nc'*rkc+rkc*nc'*rz)*rkc_norm^2 - ...
            5*(rz'*rkc)*rkc*nc'*rkc)./rkc_norm^7 +...
            ((rz'*rkc)*nc)./rkc_norm^5 );
    end
end

end


function  [centers, orientations] = projectOnCylinderSurface(positions, constraint, normCylAxis)
% projects the points at "positions" to the cylindrical constraint defined
% by "constraint" with the unit cylinder axis "normCylAxis"

    % project all InitialPositions on the cylinder constraint
    vecConstraintToCoil = bsxfun(@minus, positions, constraint.bottomCoord); % vectors from constraint origin to coil positions
    % calculate projected coil positions on the cylinder axis
    projCoilOnAxis = bsxfun(@plus, constraint.bottomCoord, normCylAxis*normCylAxis'*vecConstraintToCoil);
    % normalized vectors from projected points to coil positions
    vecAxisToCoil = normcols(positions - projCoilOnAxis);
    % map points to cylinder surface
    centers = projCoilOnAxis + constraint.radius*vecAxisToCoil;
    orientations = -vecAxisToCoil;
end

