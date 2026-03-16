function phantom = interpolateCuboidalPhantomAssembly(mrxisetup, cuboidalExtension, MNPconcentration)
% phantom = INTERPOLATECUBOIDALPHANTOMASSEMBLY(mrxisetup, cuboidalExtension, MNPconcentration)
%
% generates a discretized MNP phantom "phantom" based on the voxel grid
% defined in "mrxisetup", the cuboids defined in "cuboidalExtensions" and
% the MNP concentrations defined in "MNPconcentrations"
% ATTENTION: works only if mrxisetup.ROIData.voxelDims is available.
%
% INPUT:
% mrxisetup - standard mrxi SETUP object
% 
% cuboidalExtension - nCuboids x 1 cell; every cell holds either a 3x2
% matrix or a 3x3 matrix. The first column defines the cuboid's center of
% mass. The second column defines the x-, y-, and z-extensions of the
% cuboid. The third (optional) column defines the rotation around the x-,
% y-, and z-axes of the cuboid in rad, where the first column of the cuboid
% definition serves as the center point of the rotation.
% 
% MNPconcentration - scalar or nCuboids x 1 vector; holds the MNP
% concentrations for every cuboid defined via "cuboidalExtensions". If
% "MNPconcentrations" is a scalar, this MNP concentration is used for every
% cuboid.
% 
% OUTPUT:
% phantom - nVoxels x 1 vector; MNP phantom with interpolated MNP masses


%%% exceptions and unifying input
% check if cells are passed, otherwise convert
if ~iscell(cuboidalExtension)
    temp = cuboidalExtension;
    cuboidalExtension = cell(1);
    cuboidalExtension{1} = temp;
end

% get number of passed cuboids
nCub = length(cuboidalExtension);

% check sizes of cuboid definitions
for i = 1:nCub
    s = size(cuboidalExtension{i});
    if ~(s(1) == 3 && (s(2)==2 || s(2)==3))
        estring = ['cuboidalExtension entry ' num2str(i) ' definition wrong.'];
        error(estring);
    end
end

% check if size of passed concentrations are correctly sized
nConc = length(MNPconcentration);
if ~(nConc==1 || nConc==nCub)
    error('Wrong size of MNPconcentration,');
end
% resize MNPconcentrations if nConc=1 ~= nCub
if nConc~=nCub
    MNPconcentration = ones(nCub,1)*MNPconcentration;
end


%%% code
% get voxel data
voxelDims = mrxisetup.ROIData.voxelDims;
[~, voxRotMat] = rotatepoints([0 0 0]', mrxisetup.ROIData.rotvector);
voxExtensionsOrigin = [-voxelDims./2 voxelDims./2];
% generate voxel corner coordinates at origin
[X,Y,Z] = meshgrid( voxExtensionsOrigin(1,:)', ...
                    voxExtensionsOrigin(2,:)', ...
                    voxExtensionsOrigin(3,:)');
coordsVoxOrigin = voxRotMat * [X(:) Y(:) Z(:)]';
% make grid of sphere centers inside voxel for quick collision checks with
% cuboidal elements (only necessary if voxels are not cubic).
minDim = min(voxelDims);
nVoxCollisionSphereCenters = ceil(voxelDims./minDim);
minHalfDim = minDim/2;
voxCollisionRadius = sqrt(minHalfDim^2*3); % spherical collision check radius
% make grid of centers if necessary
if any(nVoxCollisionSphereCenters > 1)
    % make grid points in every dimension
    xDisc = linspace(   voxExtensionsOrigin(1,1)+minHalfDim, ...
                        voxExtensionsOrigin(1,2)-minHalfDim, ...
                        nVoxCollisionSphereCenters(1));
    yDisc = linspace(   voxExtensionsOrigin(2,1)+minHalfDim, ...
                        voxExtensionsOrigin(2,2)-minHalfDim, ...
                        nVoxCollisionSphereCenters(2));
    zDisc = linspace(   voxExtensionsOrigin(3,1)+minHalfDim, ...
                        voxExtensionsOrigin(3,2)-minHalfDim, ...
                        nVoxCollisionSphereCenters(3));
    % mesh grid and rotate points
    [X,Y,Z] = meshgrid( xDisc(:),yDisc(:),zDisc(:) );
    voxCollisionSphereCenters = voxRotMat*[X(:) Y(:) Z(:)]';
else
    voxCollisionSphereCenters = [0 0 0]';
end

% initialize phantom
nVox = mrxisetup.getNumberOfVoxels;
phantom = zeros(nVox,1);
% for every cuboidal element...
for i = 1:nCub
    % get current cuboid
    cub = cuboidalExtension{i};
    if size(cub,2) == 3 && ~all(cub(:,3)==0)
        rotFlag = true;
        R = getRotationMatrix(cub(1,3), cub(2,3), cub(3,3));
    else
        rotFlag = false;
    end
    
    % % make collision check of cuboidal element with voxels (like above
    % % for voxels)
    cubDims = cub(:,2);
    minDim = min(cubDims);
    nCubCollisionSphereCenters = ceil(cubDims./minDim);
    minHalfDim = minDim/2;
    cubCollisionRadius = sqrt(minHalfDim^2*3); % spherical collision check radius
    % make grid of centers if necessary
    cubExtensionsOrigin = [-cub(:,2)./2 cub(:,2)./2];
    if any(nCubCollisionSphereCenters > 1)
        % make grid points in every dimension
        xDisc = linspace(   cubExtensionsOrigin(1,1)+minHalfDim, ...
                            cubExtensionsOrigin(1,2)-minHalfDim, ...
                            nCubCollisionSphereCenters(1));
        yDisc = linspace(   cubExtensionsOrigin(2,1)+minHalfDim, ...
                            cubExtensionsOrigin(2,2)-minHalfDim, ...
                            nCubCollisionSphereCenters(2));
        zDisc = linspace(   cubExtensionsOrigin(3,1)+minHalfDim, ...
                            cubExtensionsOrigin(3,2)-minHalfDim, ...
                            nCubCollisionSphereCenters(3));
        % mesh grid and rotate points
        [X,Y,Z] = meshgrid( xDisc(:),yDisc(:),zDisc(:) );
        if rotFlag
            cubCollisionSphereCenters = R*[X(:) Y(:) Z(:)]';
        else
            cubCollisionSphereCenters = [X(:) Y(:) Z(:)]';
        end
    else
        cubCollisionSphereCenters = [0 0 0]';
    end
    
    % % check for collisions of the cuboid and the voxel spheres:
    maxCollDist = voxCollisionRadius + cubCollisionRadius;
    % add voxel collision grid points to voxel centers
    voxCheckPoints = bsxfun(@plus, mrxisetup.ROIData.voxels, ...
        permute(voxCollisionSphereCenters, [1 3 2]));
    % add cuboid collision grid points to cuboid center
    cubCheckPoints = permute(bsxfun(@plus, cub(:,1), ...
        cubCollisionSphereCenters), [1 4 3 2]);
    % calculate distances from all cuboid check points to all voxel check
    % points and check if they are shorter than maxCollDist
    if size(cubCollisionSphereCenters, 2) > 1 || size(voxCheckPoints,3) > 1
        voxUseIdxBool = any( squeeze( rssq( bsxfun(@minus, voxCheckPoints, cubCheckPoints), 1))...
                    < maxCollDist, [2 3]);
    else
        voxUseIdxBool = rssq( bsxfun(@minus, voxCheckPoints, cubCheckPoints), 1) < maxCollDist;
    end
    % if voxUseIdx is empty, continue (no collisions detected)
    if ~any(voxUseIdxBool)
        continue
    end
    voxUseIdx = 1:mrxisetup.getNumberOfVoxels;
    voxUseIdx(~voxUseIdxBool) = [];
    
    % mesh grid from cuboidal element input
    [X,Y,Z] = meshgrid( cubExtensionsOrigin(1,1:2)', ...
                        cubExtensionsOrigin(2,1:2)', ...
                        cubExtensionsOrigin(3,1:2)');
    coordsCub = [X(:) Y(:) Z(:)];
    % rotate coordinates around center of mass, if necessary
    if rotFlag
        coordsCub = (R*coordsCub')';
    end
    % add cuboid offset
    coordsCub = bsxfun(@plus, coordsCub, cub(:,1)');
    % make convex hull from cuboidal corner points
    kCub = convhull(coordsCub(:,1), coordsCub(:,2), coordsCub(:,3));
    
    % do the same as above for all voxels where a collision was detected
    for v = voxUseIdx
        % get voxel corner coordinates
        coordsVox = bsxfun(@plus, mrxisetup.ROIData.voxels(:,v), coordsVoxOrigin)';
        % make convex hull from voxel corner points
        kVox = convhull(coordsVox(:,1), coordsVox(:,2), coordsVox(:,3));
        
%         %%%%%%%%%%%% DEBUGGING:
%         figure(111); clf;
%         trisurf(kCub,coordsCub(:,1),coordsCub(:,2),coordsCub(:,3),'FaceColor','Cyan','FaceAlpha',0.2);
%         hold on;
%         trisurf(kVox,coordsVox(:,1),coordsVox(:,2),coordsVox(:,3),'FaceColor','Red','FaceAlpha',0.2);
%         daspect([1 1 1]);
%         drawnow;
%         %%%%%%%%%%%%
        
        % check if voxel corner points lie inside cuboid
        in = inhull(coordsVox,coordsCub,kCub);
        inpointsCub = coordsVox(in,:);
    
        % check if there are intersections of voxel edges with cuboid
        orig = kron(kron(coordsVox([1 4 6 7],:), ones(3,1)), ones(12,1));
        dir = kron(coordsVox([2 3 5 2 3 8 2 5 8 3 5 8], :), ones(12,1)) - orig;
        vert0 = repmat(coordsCub(kCub(:,1),:),12,1);
        vert1 = repmat(coordsCub(kCub(:,2),:),12,1);
        vert2 = repmat(coordsCub(kCub(:,3),:),12,1);
        [intersect, ~, ~, ~, xcoor] = TriangleRayIntersection (...
          orig, dir, vert0, vert1, vert2, 'lineType', 'segment', 'eps', 1e-15);
        intersectpointsCub = xcoor(intersect,:);
        
        % check if cuboid corner points lie inside voxel
        in = inhull(coordsCub,coordsVox,kVox);
        inpointsVox = coordsCub(in,:);
        
        % check if there are intersections of cuboid edges with voxel
        orig = kron(kron(coordsCub([1 4 6 7],:), ones(3,1)), ones(12,1));
        dir = kron(coordsCub([2 3 5 2 3 8 2 5 8 3 5 8], :), ones(12,1)) - orig;
        vert0 = repmat(coordsVox(kVox(:,1),:),12,1);
        vert1 = repmat(coordsVox(kVox(:,2),:),12,1);
        vert2 = repmat(coordsVox(kVox(:,3),:),12,1);
        [intersect, ~, ~, ~, xcoor] = TriangleRayIntersection (...
          orig, dir, vert0, vert1, vert2, 'lineType', 'segment', 'eps', 1e-15);
        intersectpointsVox = xcoor(intersect,:);
        
        % get unique collection of inside lying and intersection points
        % (rounded to nanometers)
        cornerpoints = unique(round([inpointsCub; intersectpointsCub; inpointsVox; intersectpointsVox], 9), 'rows');
        
%         %%%%%%%%%%%% DEBUGGING:
%         scatter3f(cornerpoints,'filled','g');
%         drawnow;
%         %%%%%%%%%%%%
        
        if ~isempty(cornerpoints) && size(cornerpoints,1)>3
            %%% check if points are collinear/coplanar
            testvecs = normdim(bsxfun(@minus, cornerpoints(2:end,:), cornerpoints(1,:)), 2);
            
            % check if points are collinear
            nTestvecs = size(testvecs,1);
            perp = cross(repmat(testvecs(1,:), nTestvecs-1, 1), testvecs(2:end,:));
            idxFirstNonzeroPerp = find(any(abs(perp)>10*eps,2), 1, 'first');
            if isempty(idxFirstNonzeroPerp)
                % the points are collinear if all cross products of the
                % normalized connection vectors between the first and all
                % other points are zero
                continue
            end
            
            % check if points are coplanar
            perp = perp(idxFirstNonzeroPerp, :);
            if all( abs(testvecs*perp')<10*eps )
                % the points are coplanar if the scalar products of all
                % connection vectors with some nonzero crossproduct are
                % zero
                continue
            end
            % calculate volume of union of voxel and cuboid. Return
            % volume=0 for coplanar/collinear points
            try
                [~, volUnion] = convhull(cornerpoints(:,1), cornerpoints(:,2), cornerpoints(:,3));
            catch ME
                if (strcmp(ME.identifier,'MATLAB:convhull:EmptyConvhull3DErrId'))
                    volUnion = 0;
                end
            end
            
%             %%%%%%%%%%%% DEBUGGING:
%             trisurf(k_corner,cornerpoints(:,1),cornerpoints(:,2),cornerpoints(:,3),...
%                 'FaceColor','Black','FaceAlpha',1);
%             xlabel x; ylabel y; zlabel z;
%             %%%%%%%%%%%%
            
            % add corresponding particle mass to respective voxel
            phantom(v) = phantom(v) + volUnion*MNPconcentration(i);
        end
    end
end






