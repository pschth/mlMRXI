classdef SETUP < handle
    % The setup class defines an MRXI setup consisting of the voxel grid
    % inside the region of interest (ROI), the excitation coil arrangement
    % as well as the sensor grid positioning.
    
%%%%%%%%%%%%%%%%%%%%%%%%%% PUBLIC MEMBERS
    properties
        
        ROIData; % contains all properties of the region of interest (ROI)
        coilData; % contains all coil properties
        sensorData; % constains all sensor properties
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%% PUBLIC FUNCTIONS
    methods (Access = public)
        
        % constructor
        function o = SETUP()
            % constructor
        end
        function setupcopy = copySetup(o)
            % copies the entire settings into the new object "setupcopy".
            % (By copying an object via equal sign, changes conducted in
            % copied object also occur in original object. "copySetup"
            % prevents that.)
            setupcopy = SETUP();
            setupcopy.ROIData = o.ROIData;
            setupcopy.coilData = o.coilData;
            setupcopy.sensorData = o.sensorData;
        end
        % setter functions
        function o = setVoxelGrid(o, ROI, resolution, offset, rotvector, shape)
            % constructs an equally spaced voxel grid with the resolution
            % "resolution" in the boundary region "ROI"
            %
            % INPUT: ROI - 3x2 matrix, contains the boundary coordinates in
            % x and y and z direction (i.e. [-x +x; -y +y; -z +z])
            %
            % resolution - 3x1 vector, contains the number of voxels in x
            % and y and z direction
            %
            % offset - 3x1 vector (optional), ROI offset from origin;
            % Default: zeros(3,1)
            %
            % rotvector - 3x1 vector (optional), ROI rotation direction;
            % Default: [0 0 1].
            %
            % shape - string (optional), defines shape of the voxel grid.
            % Can be 'cuboid', 'cylinder' or 'sphere'. The radius/radii of
            % the cylinder is defined by the dimensions of the ROI in x-
            % and y-direction. The radius/radii of the sphere is defined by
            % the dimensions of the ROI in x-, y-, and z-direction.
            % Default: 'cuboid'
            
            %%% exceptions
            if ~(numel(ROI) == 6 && numel(resolution) == 3)
                error('Wrong dimension of input variables.');
            end
            
            % transpose ROI if necessary
            if size(ROI,1) == 2 && size(ROI,2) == 3
                ROI = ROI';
            end

            if ~all(size(ROI)==[3 2])
               error('ROI has the wrong size.');
            end

            if nargin < 4 || isempty(offset)
                offset = zeros(3,1);
            end
            offset = offset(:);
            if ~all(size(offset)==[3 1])
                error('offset has the wrong size.');
            end
            
            if nargin < 5
                rotvector = [0 0 1];
            elseif isempty(rotvector)
                rotvector = [0 0 1];
            end
            rotvector = normcols(rotvector(:));
            if ~all(size(rotvector)==[3 1])
                error('rotvector has the wrong size.');
            end
            
            if nargin < 6
                shape = 'cuboid';
            elseif isempty(shape)
                shape = 'cuboid';
            end
            if ~(strcmp(shape,'cuboid') || ...
                 strcmp(shape,'cylinder') || ...
                 strcmp(shape,'sphere')) 
                error('Shape not defined.');
            end
            
            
            %%% code
            % get initial offset of ROI and remove it from ROI to enable
            % clean transformations
            initOffset = mean(ROI, 2);
            ROI = bsxfun(@minus, ROI, initOffset);
            
            % calculate voxel coordinates based on ROI and resolution
            x = linspace(ROI(1,1),ROI(1,2),resolution(1)*2+1);
            x = x(2:2:end-1);
            y = linspace(ROI(2,1),ROI(2,2),resolution(2)*2+1);
            y = y(2:2:end-1);
            z = linspace(ROI(3,1),ROI(3,2),resolution(3)*2+1);
            z = z(2:2:end-1);
            
            % make volumetric grid from voxel coordinates
            [Y,X,Z] = meshgrid(y,x,z);
            
            % for phantom visualization; defines which voxels are marked
            % for deletion. (the phantom visualization grid is always
            % defined as cuboidal grid. delIdx defines voxels that shall
            % not be displayed in case of for instance cylindrical or
            % spherical shape of the ROI.)
            delIdx = false(o.getNumberOfVoxels,1);
            
            % delete voxels outside of cylindrical shape if defined
            if strcmp(shape,'cylinder')
                % if shape is cylinder, remove all voxels outside of cylinder
                % radius along z-axis
                Xc = X(:) - mean(X(:)); % centered x-coordinates
                Yc = Y(:) - mean(Y(:)); % centered y-coordinates
                rx = diff(ROI(1,:))/2; % radius in x-direction
                ry = diff(ROI(2,:))/2; % radius in y-direction
                Yc = (rx/ry) .* Yc; % stretch y-coordinates, so that radii in x- and y-direction are equal (=rx)
                R = sqrt(Xc.^2 + Yc.^2); % calculate radius to central axis at every point
                delIdx = R>rx; % mark every coordinate greater rx for deletion
                X = X(~delIdx); % delete voxels outside of shape
                Y = Y(~delIdx); % delete voxels outside of shape
                Z = Z(~delIdx); % delete voxels outside of shape
            end
            % delete voxels outside of spherical shape if defined
            if strcmp(shape,'sphere')
                % if shape is sphere, remove all voxels outside of radius
                Xc = X(:) - mean(X(:)); % centered x-coordinates
                Yc = Y(:) - mean(Y(:)); % centered y-coordinates
                Zc = Z(:) - mean(Z(:)); % centered z-coordinates
                rx = diff(ROI(1,:))/2; % radius in x-direction
                ry = diff(ROI(2,:))/2; % radius in y-direction
                rz = diff(ROI(3,:))/2; % radius in z-direction
                Yc = (rx/ry) .* Yc; % stretch y-coordinates, so that radii in x- and y-direction are equal (=rx)
                Zc = (rx/rz) .* Zc; % stretch Z-coordinates, so that radii in x- and z-direction are equal (=rx)
                R = sqrt(Xc.^2 + Yc.^2 + Zc.^2); % calculate radius to central axis at every point
                delIdx = R>rx; % mark every coordinate greater rx for deletion
                X = X(~delIdx); % delete voxels outside of shape
                Y = Y(~delIdx); % delete voxels outside of shape
                Z = Z(~delIdx); % delete voxels outside of shape
            end
            
            % get number of voxels
            nVoxel = numel(X);
            
            % define voxel grid as 3 x nVoxel coordinate matrix
            o.ROIData.voxels = [X(:) Y(:) Z(:)]';
            
            % rotate voxel grid from [0 0 1] to new orientation 'rotvector'
            o.ROIData.voxels = rotatepoints(o.ROIData.voxels, rotvector);
            
            % translate voxel grid back to initial ROI offset plus the
            % additional input parameter 'offset'
            o.ROIData.voxels = o.ROIData.voxels + repmat(initOffset+offset,1,nVoxel);
            
            % % % set remaining ROIData
            o.ROIData.voxelDims = diff(ROI,[],2)./resolution(:); % calculate size of one voxel
            o.ROIData.ROI = ROI; % set ROI without (initial) offset
            o.ROIData.resolution = resolution; % set ROI resolution
            o.ROIData.offset = initOffset+offset; % set offset property as initial ROI offset plus input parameter offset
            o.ROIData.rotvector = rotvector; % set ROI orientation
            o.ROIData.shape = shape; % set ROI shape
            o.ROIData.voxDelIdx = delIdx; % store voxel indices marked for deletion
            
        end
        function o = setVoxelCoordinates(o, voxels)
            % the voxel coordinates can be set manually. The ROI is
            % calculated automatically. (Automatic detection of remaining
            % ROI parameters works only well for voxel grids regularly
            % spaced along the x-/y-/z-axes.
            %  
            % INPUT:
            % voxels - 3 x nVoxel matrix; centerpoints of the coils
            
            %%% exceptions
            % check if matrices have right sizes
            if ~any(size(voxels) == 3)
                error('Wrong size of "voxels" matrix.');
            end
            
            % transpose matrices if necessary
            if size(voxels,1) ~= 3 && size(voxels,2) == 3
                voxels = voxels';
            end
            
            %%% code
            % set voxel coordinates
            o.ROIData.voxels = voxels;
            
            % estimate resolution from voxel coordinates
            o.getResolutionFromVoxels(voxels);

            % get unique voxel positions in x, y and z
            x = unique(voxels(1,:));
            y = unique(voxels(2,:));
            z = unique(voxels(3,:));
            
            % estimate voxel dimensions
            if length(x) == 1
                dx = 0;
            else
                dx = x(2)-x(1);
            end
            if length(y) == 1
                dy = 0;
            else
                dy = y(2)-y(1);
            end
            if length(z) == 1
                dz = 0;
            else
                dz = z(2)-z(1);
            end
            
            % calculate voxel distances in x, y and z
            centerdist = [dx; dy; dz];
            
            % calculate ROI
            o.ROIData.ROI = [   min(x)-centerdist(1)/2 max(x)+centerdist(1)/2;...
                                min(y)-centerdist(2)/2 max(y)+centerdist(2)/2;...
                                min(z)-centerdist(3)/2 max(z)+centerdist(3)/2];
            % calculate offset
            o.ROIData.offset = mean(o.ROIData.ROI,2);
            
            % correct ROI with offset
            o.ROIData.ROI = [o.ROIData.ROI(:,1) - o.ROIData.offset ...
                             o.ROIData.ROI(:,2) - o.ROIData.offset];
                         
            % set rotation vector
            o.ROIData.rotvector = [0 0 1]';
            
            % remove voxel dimensions if necessary (for stable
            % setup visualization)
            if isfield(o.ROIData, 'voxelDims')
                o.ROIData = rmfield(o.ROIData, 'voxelDims');
            end
            
        end
        function o = setNewResolution(o, resolution)
            % constructs an equally spaced voxelgrid with the dimension
            % "resolution" in the boundary region ROI. ROI has to be 
            % defined for this function to work. 
            %
            % INPUT:
            % resolution - 3x1 vector, contains the number of voxels
            % in x and y and z direction
            
            %%% exceptions
            if ~(numel(resolution) == 3)
                error('Wrong dimension of input variables.');
            end
            if isempty(o.ROIData)
                error('Undefined ROI.')
            elseif ~isfield(o.ROIData,'ROI')
                error('Undefined ROI.')
            end
            
            %%% code
            % set voxel grids using new resolution (also supports older
            % versions of the code)
            % Info for pyMRXI: only the first 'if' entry is needed for
            % translating the code to python
            if isfield(o.ROIData, 'shape')
                o.setVoxelGrid(o.ROIData.ROI, resolution, o.ROIData.offset, o.ROIData.rotvector, o.ROIData.shape);
            elseif isfield(o.ROIData, 'offset')
                o.setVoxelGrid(o.ROIData.ROI, resolution, o.ROIData.offset, o.ROIData.rotvector);
            else
                o.setVoxelGrid(o.ROIData.ROI, resolution);
            end
        end
        function o = setCoils(o, centers, orientations, coilpattern, coilpatternassignment)
            % specifies the excitation coil positions, orientations, and
            % shapes
            %
            % INPUT:
            % centers - 3 x nCoils matrix; centerpoints of the coils
            %
            % orientations - 3 x nCoils matrix; normalvectors of the coils
            %
            % coilpattern - (optional) 3 x nCoils matrix OR 1 x
            % nCoilPatterns cell containing nCoilPatterns 3 x nCoils
            % matrices;
            % specifies the points approximating the coil design. If not 
            % stated, the coils will be defined as magnetic dipoles. Coil
            % patterns are assumed to be centered at origin and oriented in
            % z-direction.
            %
            % coilpatternassignment - (optional) nCoils vector; contains
            % for each coil the coilpattern index for the coilpattern they
            % are assigned to. Must be specified for the use of
            % nCoilPatterns>1 coil patterns. Default: ones(nCoils,1)
            
            % empty old coilData except for coil constraints which is not
            % specified by 'setCoils'
            if ~isempty(o.coilData)
                coilDataFieldNames = fieldnames(o.coilData);
                for fn = 1:length(coilDataFieldNames)
                    if ~strcmp(coilDataFieldNames{fn},'constraints')
                        o.coilData = rmfield(o.coilData, coilDataFieldNames{fn});
                    end
                end
            end
            
            %%% exceptions
            % check if matrices have right sizes
            if ~any(size(centers) == 3)
                error('Wrong size of "centers" matrix.');
            end
            if ~any(size(orientations) == 3)
                error('Wrong size of "orientations" matrix.');
            end
            
            % transpose matrices if necessary
            if size(centers,1) ~= 3 && size(centers,2) == 3
                centers = centers';
            end
            if size(orientations,1) ~= 3 && size(orientations,2) == 3
                orientations = orientations';
            end
            
            % check if coils are dipoles or have patterns
            if nargin < 4
                FLAG_isDipole = true;
            elseif isempty(coilpattern)
                FLAG_isDipole = true;
            else
                % coil is no dipole -> check matrix size and transpose if
                % necessary
                if iscell(coilpattern)
                    % if multiple coil patterns are defined, check all coil
                    % patterns
                    for i = 1:numel(coilpattern)
                        if ~any(size(coilpattern{i}) == 3)
                            errstring = ['Wrong size of "coilpattern" matrix ' num2str(i) '.'];
                            error(errstring);
                        end
                        if size(coilpattern{i},1) ~= 3 && size(coilpattern{i},2) == 3
                            coilpattern{i} = coilpattern{i}';
                        end
                    end
                else
                    % if one coil pattern is defined, check it
                    if ~any(size(coilpattern) == 3)
                        error('Wrong size of "coilpattern" matrix.');
                    end
                    if size(coilpattern,1) ~= 3 && size(coilpattern,2) == 3
                        coilpattern = coilpattern';
                    end
                    % insert single coil pattern into cell to guarantee
                    % uniform input for subsequent computations
                    coilpattern_temp = coilpattern;
                    clear coilpattern;
                    coilpattern{1} = coilpattern_temp;
                end
                FLAG_isDipole = false; % defined coils are no dipoles
            end
            
            % check if coilpatternassignment was stated and SETUP does not
            % consist of dipolar excitation sources
            if nargin < 5 && ~FLAG_isDipole
                coilpatternassignment = ones(size(centers,2),1);
            elseif ~FLAG_isDipole && isempty(coilpatternassignment)
                coilpatternassignment = ones(size(centers,2),1);
            elseif ~FLAG_isDipole
                if ~numel(coilpatternassignment) == size(centers,2)
                    error('Wrong size of "coilpatternassignment" vector.');
                end
            end
            
            %%% code
            % set object properties to input data
            o.coilData.centers = centers;
            orientations = normcols(orientations); % normalize orientations to unit vector lengths
            o.coilData.orientations = orientations;
            o.coilData.isDipole = FLAG_isDipole;
            
            % if coil is no dipole, rotate and translate the coil pattern
            % according to centers and orientations of coils
            if ~FLAG_isDipole
                nCoils = size(centers,2); % get number of coils
                % fill coilData of SETUP object
                o.coilData.coilpattern = coilpattern;
                o.coilData.coilpatternassignment = coilpatternassignment;
                % initialize rotation matrices of the coils
                o.coilData.rotMatrices = zeros(3,3,nCoils);
                
                % for each coil, translate and rotate the coil pattern
                % according to 'centers' and 'orientations' and store the
                % final locations of the coil pattern's filamentary
                % segments under coilData.segments.c[coil index]
                for c = 1:nCoils
                    % generate structure name according to coil index
                    structname = ['c' num2str(c)];
                    % rotate assigned coil pattern to new orientation
                    [rotatedcoilpattern, o.coilData.rotMatrices(:,:,c)] = ...
                        rotatepoints(coilpattern{coilpatternassignment(c)}, orientations(:,c));
                    % translate the assigned rotated coil pattern to new
                    % center
                    o.coilData.segments.(structname) = ...
                        bsxfun(@plus, centers(:,c), rotatedcoilpattern);                    
                end
            end
        end
        function o = addCoils(o, centers, orientations, coilpattern, coilpatternassignment)
            % adds excitation coils. Can only be used if "setCoils" was
            % called at least once before.
            %
            % INPUT:
            % centers - 3 x nCoils matrix; centerpoints of the coils
            %
            % orientations - 3 x nCoils matrix; normalvectors of the coils
            %
            % coilpattern - (optional) 3 x nCoils matrix OR 1 x
            % nCoilPatterns cell containing nCoilPatterns 3 x nCoils
            % matrices;
            % specifies the points approximating the coil design. If not 
            % stated, the coils will be defined as magnetic dipoles. Coil
            % patterns are assumed to be centered at origin and oriented in
            % z-direction.
            %
            % coilpatternassignment - (optional) nCoils vector; contains
            % for each coil the coilpattern index for the coilpattern they
            % are assigned to. Must be specified for the use of
            % nCoilPatterns>1 coil patterns. Default: ones(nCoils,1)
            
            %%% exceptions  
            % check if setCoils was called at least once before
            if ~isfield(o.coilData, 'centers') || ~isfield(o.coilData, 'orientations')
                error('Use the function "setCoils" initially before adding coils with "addCoils".')
            end
            
            % check if matrices have right sizes
            if ~any(size(centers) == 3)
                error('Wrong size of "centers" matrix.');
            end
            if ~any(size(orientations) == 3)
                error('Wrong size of "orientations" matrix.');
            end
            
            % transpose matrices if necessary
            if size(centers,1) ~= 3 && size(centers,2) == 3
                centers = centers';
            end
            if size(orientations,1) ~= 3 && size(orientations,2) == 3
                orientations = orientations';
            end
            
            % check if coils are dipoles or have patterns and the initial
            % setup has the same properties
            if nargin > 3
                if o.coilData.isDipole
                    error('Initial SETUP coils were specified as dipoles. Cannot be paired with non-dipoles.')
                end
                % check matrix size(s) and transpose if necessary
                if isempty(coilpattern)
                    coilpattern = [];
                elseif iscell(coilpattern)
                    for i = 1:numel(coilpattern)
                        if ~any(size(coilpattern{i}) == 3)
                            errstring = ['Wrong size of "coilpattern" matrix ' num2str(i) '.'];
                            error(errstring);
                        end
                        if size(coilpattern{i},1) ~= 3 && size(coilpattern{i},2) == 3
                            coilpattern{i} = coilpattern{i}';
                        end
                    end
                else
                    if ~any(size(coilpattern) == 3)
                        error('Wrong size of "coilpattern" matrix.');
                    end
                    if size(coilpattern,1) ~= 3 && size(coilpattern,2) == 3
                        coilpattern = coilpattern';
                    end
                    % insert single coil pattern into cell to guarantee
                    % uniform input for subsequent computations
                    coilpattern_temp = coilpattern;
                    clear coilpattern;
                    coilpattern{1} = coilpattern_temp;
                end
            end
            
            % check if coilpatternassignment was stated
            if nargin < 5 && ~o.coilData.isDipole
                coilpatternassignment = ones(size(centers,2),1);
            elseif isempty(coilpatternassignment) && ~isDipole
                coilpatternassignment = ones(size(centers,2),1);
            elseif ~o.coilData.isDipole
                if ~any(size(coilpatternassignment) == size(centers,2))
                    error('Wrong size of "coilpatternassignment" vector.');
                end
            end
            
            %%% code
            % append input data to object properties
            nCoilsOld = size(o.coilData.centers,2);
            o.coilData.centers = [o.coilData.centers centers];
            o.coilData.orientations = [o.coilData.orientations normcols(orientations)];
            
            % if coil is no dipole, rotate and translate the coil pattern
            % according to centers and orientations of coils
            if ~o.coilData.isDipole
                % get coil indices of added coils
                nCoilsNew = size(centers,2);
                newCoilIdx = nCoilsOld+1:nCoilsOld+nCoilsNew;
                % append new coil patterns and coilpatternassignments
                o.coilData.coilpattern = [o.coilData.coilpattern(:); coilpattern(:)];
                o.coilData.coilpatternassignment = [o.coilData.coilpatternassignment(:); coilpatternassignment(:)];
                % initialize new rotation matrices
                o.coilData.rotMatrices(:,:,newCoilIdx) = zeros(3,3,nCoilsNew);
                
                % for each coil...
                for c = newCoilIdx
                    % generate structure name according to coil index
                    structname = ['c' num2str(c)];
                    % rotate assigned coil pattern to new orientation
                    [rotatedcoilpattern, o.coilData.rotMatrices(:,:,c)] = ...
                        rotatepoints(o.coilData.coilpattern{o.coilData.coilpatternassignment(c)}, o.coilData.orientations(:,c));
                    % translate the assigned rotated coil pattern to new
                    % center
                    o.coilData.segments.(structname) = ...
                        bsxfun(@plus, o.coilData.centers(:,c), rotatedcoilpattern);                    
                end  
            end
        end
        function o = setCoilCenters(o, newCenters)
            % moves existing coils to new positions
            %
            % INPUT:
            % centers - 3 x nCoils matrix; centerpoints of
            % the coils
            
            %%% exceptions
            % check if matrices have right sizes
            if ~any(size(newCenters) == 3)
                error('Wrong size of "centers" matrix.');
            end
            
            % transpose matrices if necessary
            if size(newCenters,1) ~= 3 && size(newCenters,2) == 3
                newCenters = newCenters';
            end
            
            %%% code
            % reset segments if necessary
            if ~o.coilData.isDipole
                for c = 1:o.getNumberOfCoils
                    % get structure name according to coil index
                    coilName = ['c' num2str(c)];
                    % move filamentary segments to new position
                    o.coilData.segments.(coilName) = bsxfun(@plus, ...
                        o.coilData.segments.(coilName), ...
                        newCenters(:,c) - o.coilData.centers(:,c));
                end                
            end
            % set coil centers
            o.coilData.centers = newCenters;
        end
        function o = setCoilRadii(o, radii)
            % scales coil radii to given radii-vector. Should only be used
            % for cylindrical coils. Can be called after setCoils function
            % with non-dipole coils has been called.
            %
            % INPUT:
            % radii - nCoils vector; absolute values of cylindrical coil
            % radii
           
            %%% exceptions
            % check if matrices have right sizes
            if ~any(size(radii) == 1) && ~any(size(radii) == o.getNumberOfCoils)
                error('Wrong size of "radii" vector.');
            end
            
            %%% code
            % get original radius from coil pattern
            origRadius = o.getCoilPatternRadius;
            
            % for each coil...
            for c = 1:o.getNumberOfCoils
                % get structure name according to coil index
                structname = ['c' num2str(c)];
                % get relative scaling factors for coils
                scalingFactor = radii(c)/origRadius(o.coilData.coilpatternassignment(c));
                % scale, rotate and translate coil pattern to original
                % position
                o.coilData.segments.(structname) = ...
                    bsxfun(@plus,...
                    scalingFactor*o.coilData.rotMatrices(:,:,c)*o.coilData.coilpattern{o.coilData.coilpatternassignment(c)},...
                    o.coilData.centers(:,c));                 
            end
            
            
        end
        function o = setCoilCenterConstraint(o, type, varargin)
            % sets constraint for coil positions; see description of
            % 'setCenterConstraint' for more information
            
            o.setCenterConstraint('coilData', type, varargin);
        end
        function o = removeCoils(o, idx)
            % deletes coil with index 'idx' (either logical or decimal).
            %
            % INPUT:
            % idx - scalar or vector; contains index or indices (either
            % logical or decimal) of coils to delete. If idx is left blank,
            % the function deletes all coils.
            
            % determine coils to keep
            % if idx is empty, delete all coils
            if nargin < 2
                keepIdx = [];
            elseif isempty(idx)
                keepIdx = [];
            else
            % otherwise, determine coil indices to keep
            if islogical(idx)
                keepIdx = find(~idx);
            else
                keepIdx = setdiff(1:o.getNumberOfCoils, idx);
            end
            end
            
            if o.coilData.isDipole
                % keep remaining coils
                o.coilData.centers = o.coilData.centers(:, keepIdx);
                o.coilData.orientations = o.coilData.orientations(:, keepIdx);
            else
                % keep remaining coils
                o.coilData.centers = o.coilData.centers(:, keepIdx);
                o.coilData.orientations = o.coilData.orientations(:, keepIdx);
                o.coilData.rotMatrices = o.coilData.rotMatrices(:, :, keepIdx);
                o.coilData.coilpatternassignment = o.coilData.coilpatternassignment(keepIdx);
                % reorder the filamentary segment struct of the coil data
                temp = struct2cell(o.coilData.segments);
                o.coilData.segments = struct;
                for c = 1:o.getNumberOfCoils
                    cname = ['c' num2str(c)];
                    o.coilData.segments.(cname) = temp{keepIdx(c)};
                end
                % remove coil patterns and reorder coilpatternassignment
                % indices if a coil pattern is not used anymore
                [remaining_coilpatterns, ~, o.coilData.coilpatternassignment] = ...
                    unique(o.coilData.coilpatternassignment);
                o.coilData.coilpattern = o.coilData.coilpattern(remaining_coilpatterns);
            end         
        end
        function o = removeCoilCenterConstraints(o, idx)
            % deletes coil center constraints with decimal indices 'idx'.
            % If called empty, all constraints are deleted.
            
            if isfield(o.coilData, 'constraints')
                if nargin < 2
                    % remove all constraints
                    o.coilData = rmfield(o.coilData, 'constraints');
                else
                    % remove constraints with indices idx
                    if max(idx) > length(o.coilData.constraints)
                        error('Cannot delete coil constraints: specified idx too large.');
                    else
                        o.coilData.constraints(idx) = [];
                    end
                end
            end
        end
        function o = setSensors(o, centers, orientations)
            % specifies the sensor positions and the orientations of their
            % sensitive axes
            %
            % INPUT
            % centers - 3 x nCoils matrix; centerpoints of the sensors
            %
            % orientations - 3 x nCoils matrix; normalvectors of the
            % sensors
            
            %%% exceptions
            % check if matrices have right sizes
            if ~any(size(centers) == 3)
                error('Wrong size of "centers" matrix.');
            end
            if ~any(size(orientations) == 3)
                error('Wrong size of "orientations" matrix.');
            end
            
            % transpose matrices if necessary
            if size(centers,1) ~= 3 && size(centers,2) == 3
                centers = centers';
            end
            if size(orientations,1) ~= 3 && size(orientations,2) == 3
                orientations = orientations';
            end
            
            %%% code
            o.sensorData.centers = centers;
            o.sensorData.orientations = normcols(orientations); % normalize orientation vectors to unit lengths
        end
        function o = setSensorCenterConstraint(o, type, varargin)
            % sets constraint for sensor positions; see description of
            % 'setCenterConstraint' for more information
            
            o.setCenterConstraint('sensorData', type, varargin);
       end
        

 
        % getter functions
        function res = getResolution(o)
            % returns the resolution of the object
            res = o.ROIData.resolution;
        end
        function nVox = getNumberOfVoxels(o)
            % returns the number of voxels in the setup grid
            if isempty(o.ROIData)
                nVox = 0;
            elseif isempty(o.ROIData.voxels)
                nVox = 0;
            else
                nVox = size(o.ROIData.voxels,2);
            end
        end
        function nCoils = getNumberOfCoils(o)
            % returns the number of setup coils
            if isprop(o,'coilData') && isempty(o.coilData)
                nCoils = 0;
            elseif isfield(o.coilData,'centers') && isempty(o.coilData.centers)
                nCoils = 0;
            else
                nCoils = size(o.coilData.centers,2);
            end
        end
        function nSensors = getNumberOfSensors(o)
            % returns the number of setup sensors
            if isprop(o,'sensorData') && isempty(o.sensorData)
                nSensors = 0;
            elseif isfield(o.sensorData,'centers') && isempty(o.sensorData.centers)
                nSensors = 0;
            else
                nSensors = size(o.sensorData.centers,2);
            end
        end
        function radii = getCoilRadii(o)
            % returns the maximum extension of the coils along the plane
            % perpendicular to the orientation vector of the coil
            if o.coilData.isDipole
                error('Coils are defined as dipoles. Cannot retrieve coil radii.');
            end
            
            % initialize radius vector
            nCoils = o.getNumberOfCoils;
            radii = zeros(nCoils,1);

            for c = 1:nCoils
                % get structure name of filamentary segments
                structname = ['c' num2str(c)];
                % get distance vectors of filament segments to center of
                % coil
                distVecs = bsxfun(@minus, ...
                    o.coilData.segments.(structname), o.coilData.centers(:,c));
                % get distances of filament segments to center of coil
                % (hypotenuse)
                hyps = rssq(distVecs);
                % get projections of distVecs on orientation vector
                % (adjacent leg)
                adjs =  o.coilData.orientations(:,c)' * distVecs;
                % calculate distances between filamentary segments and
                % centers projected on the plane perpendicular to the
                % orientation
                dists = hyps .* sin(acos(adjs./hyps));
                % get maximum value of these distances (=coil radius)
                radii(c) =  max(dists);
            end
        end
        function radii = getCoilPatternRadius(o)
            % returns the maximum coil pattern extensions (=radii) measured
            % from the origin
            if o.coilData.isDipole
                error('Coils are defined as dipoles. Cannot retrieve coil radii.');
            end

            % initialize radius vector
            radii = zeros(numel(o.coilData.coilpattern),1);
            % radii are the maximum x-y-distances, since coil patterns are
            % oriented along the z-axis and centered in the origin
            for cp = 1:numel(o.coilData.coilpattern)
                radii(cp) = max(rssq(o.coilData.coilpattern{cp}(1:2,:),1));
            end
        end
        function unitRadiusCoilPattern = getCoilPatternWithUnitRadius(o)
            % returns cell with the coil patterns with unit radius and
            % oriented along z-axis
            if o.coilData.isDipole
                error('Coils are defined as dipoles. Cannot retrieve unit radius coil pattern.');
            end
            
            % get number of different coil patterns
            nCoilPatterns = numel(o.coilData.coilpattern);
            % get coil pattern radii
            cpRadii = o.getCoilPatternRadius();
            % scale patterns to unit radii
            unitRadiusCoilPattern = cell(nCoilPatterns,1);
            for cp = 1:nCoilPatterns
                unitRadiusCoilPattern{cp} = o.coilData.coilpattern{cp}./cpRadii(cp);
            end            
        end
        function positions = getCylindricalConstraintPositions(o)
            % returns the positions required for cylindrical coil placement
            % optimization
            %
            % OUTPUT
            % positions - 2 x nCoils matrix; (1) angles on and (2)
            % distances along the cylindrical constraint
            
            % check if constraint is defined properly
            if ~isfield(o.coilData, 'constraints')
                error('No constraints defined on this object.');
            elseif length(o.coilData.constraints) ~= 1
                error('Only 1 cylindrical constraint allowed.');
            elseif ~strcmp(o.coilData.constraints{1}.type, 'cylindrical')
                error('Only 1 cylindrical constraint allowed.');
            end
            constraint = o.coilData.constraints{1};
            
            % remove offset of constraint
            centers = bsxfun(@minus, o.coilData.centers, constraint.bottomCoord);
            % rotate constraint back to [0 0 1]'
            cAxis = constraint.topCoord - constraint.bottomCoord;
            centers = rotatepoints(centers, [0 0 1]', cAxis);
            % calculate angles and set axis heights
            positions(1,:) = atan(centers(2,:)./centers(1,:)) ... calculate angle between -pi/2 and + pi/2
                + (centers(1,:)<0)*pi; % and add pi if center lies in 2nd or 3rd quadrant
            positions(2,:) = centers(3,:); % axis height is just the transformed z-distance
        end
        function xyExtensions = getPlanarCoilExtensions(o)
            % returns coil extensions perpendicular to their normal axes
            
            % get informations about which and how many coil patterns were
            % used
            usedCoilPatterns = unique(o.coilData.coilpatternassignment);
            nUsedCoilPatterns = length(usedCoilPatterns);
            
            % calculate the extensions of the used coil pattern
            xyExtensionsPattern = zeros(2,nUsedCoilPatterns);
            for cp = 1:nUsedCoilPatterns
                xyExtensionsPattern(:,cp) = max(o.coilData.coilpattern{cp}(1:2,:),[],2) - ...
                                            min(o.coilData.coilpattern{cp}(1:2,:),[],2);
            end
            
            % get radii of the original coil patterns
            coilPatternRadii = o.getCoilPatternRadius;
            
            % get radii of the used coils
            coilRadii = o.getCoilRadii;
            
            % for each coil, scale the respective coil pattern extensions
            % with the ratio of coilRadii to coilPatternRadii
            nCoils = o.getNumberOfCoils;
            xyExtensions = zeros(2,nCoils);
            for c = 1:nCoils
                cp = usedCoilPatterns == o.coilData.coilpatternassignment(c);
                xyExtensions(:,c) = ...
                    xyExtensionsPattern(:,cp) * coilRadii(c)/coilPatternRadii(cp);
            end
        end
        function [roi, cornerpoints] = getTransformedROIBoundary(o)
            % returns the ROI after translation and reorientation according
            % to ROIData.offset and ROIData.rotvector
            
            % get corner points of initial cuboidal ROI
            roi = o.ROIData.ROI;
            cornerpoints = [roi(1,1) roi(1,2) roi(1,1) roi(1,2) roi(1,1) roi(1,2) roi(1,1) roi(1,2); ...
                            roi(2,1) roi(2,1) roi(2,2) roi(2,2) roi(2,1) roi(2,1) roi(2,2) roi(2,2); ...
                            roi(3,1) roi(3,1) roi(3,1) roi(3,1) roi(3,2) roi(3,2) roi(3,2) roi(3,2)];
            if isfield(o.ROIData, 'rotvector') && isfield(o.ROIData, 'offset')
                % get center of ROI (should be centered already but just to
                % be sure...)
                roicenter = mean(roi,2);

                % translate cornerpoints to origin
                cornerpoints = bsxfun(@minus, cornerpoints, roicenter);

                % rotate cornerpoints
                cornerpoints = rotatepoints(cornerpoints, o.ROIData.rotvector);

                % translate cornerpoints to final position
                cornerpoints = bsxfun(@plus, cornerpoints, roicenter+o.ROIData.offset);

                % get minimum/maximum extent of ROI
                roi = [min(cornerpoints,[],2), max(cornerpoints,[],2)];
            end
        end

        

        % setup visualization
        function [voxHandle, coilHandle, sensHandle] = visualize(o, VoxAsScatter)
            % visualizes current setup config. 
            % ATTENTION: voxels can only be displayed as cuboids if
            % ROIData.voxelDims is given. This is not automatically defined
            % if only voxel coordinates are passed (with
            % setVoxelCoordinates). If ROIData.voxelDims is not given,
            % voxels positions are represented as scatter plot.
            % 
            % INPUT:
            % VoxAsScatter (optional) - boolean; if true, voxels are
            % represented as scatter plot rather than the ROI as box;
            % default: false
            
            %%% exceptions
            if nargin < 2
                VoxAsScatter = false;
            elseif isempty(VoxAsScatter)
                VoxAsScatter = false;
            end
            
            %%% code
            [voxHandle, coilHandle, sensHandle] = visualizeSetup(o, VoxAsScatter);
            xlabel('x'); ylabel('y'); zlabel('z');
            daspect([1 1 1]);
            
            % disable axis limits
            set(gca,'XLim',[-Inf Inf]);
            set(gca,'YLim',[-Inf Inf]);
            set(gca,'ZLim',[-Inf Inf]);
        end
        function [constraintHandle] = visualizeConstraints(o)
            % visualizes the coil and sensor placement constraints
            
            %%% exceptions
            if ~isfield(o.coilData, 'constraints') && ~isfield(o.sensorData, 'constraints')
                warning('Object has no defined constraints.');
            else

                %%% code
                % initialize
                holdStatus = ishold; % store current hold status of figure
                % get number of coil constraints if there are any
                if isfield(o.coilData, 'constraints')
                    nCConstraints = length(o.coilData.constraints);
                else
                    nCConstraints = 0;
                end
                % get number of sensor constraints if there are any
                if isfield(o.sensorData, 'constraints')
                    nSConstraints = length(o.sensorData.constraints);
                else
                    nSConstraints = 0;
                end
                % initialize figure handle that will contain all constraint
                % visualizations
                constraintHandle = cell(nCConstraints+nSConstraints,1);

                % plot all constraints
                for c_or_s = 1:2 % coil constraints (=1) or sensor constraints (=2)
                    switch c_or_s
                        case 1
                            % coil constraints
                            if nCConstraints == 0
                                continue
                            end
                            color = [255, 150, 0]./255; % define color for coil constraints
                            iterEnd = nCConstraints; % get number of coil constraints for subsequent for loop
                            constraints = o.coilData.constraints; % copy coil constraints
                            fillIdxOffset = 0; % for filling 'constraintHandle'
                        case 2
                            % sensor constraints
                            if nSConstraints == 0
                                continue
                            end
                            color = [0, 105, 255]./255; % define color for sensor constraints
                            iterEnd = nSConstraints; % get number of sensor constraints for subsequent for loop
                            constraints = o.sensorData.constraints; % copy sensor constraints
                            fillIdxOffset = nCConstraints; % for filling 'constraintHandle'
                    end

                    % loop through all coil/sensor constraints
                    for i = 1:iterEnd
                        % plot constraints
                        hold on;
                        switch constraints{i}.type % 'planar' or 'cylindrical' constraints can be defined
                            case 'cylindrical'
                            % make cylinder with constraint radius and length 1
                            [X,Y,Z] = cylinder(constraints{i}.radius,50);
                            % stretch cylinder to correct length
                            cAxis = constraints{i}.topCoord - constraints{i}.bottomCoord; % central axis vector of cylinder
                            lAxis = norm(cAxis); % length of central axis vector
                            ncAxis = cAxis./lAxis; % central vector with length = 1
                            Z = Z.*lAxis; % scale cylinder z-coordinates to correct length
                            % rotate to constraint axis and translate to constraint
                            % origin
                            XYZ = rotatepoints([X(:) Y(:) Z(:)]', ncAxis); % rotation
                            X = reshape(XYZ(1,:), size(X)) + constraints{i}.bottomCoord(1); % translation x
                            Y = reshape(XYZ(2,:), size(X)) + constraints{i}.bottomCoord(2); % translation y
                            Z = reshape(XYZ(3,:), size(X)) + constraints{i}.bottomCoord(3); % translation z
                            % plot cylinder and store surf handle
                            constraintHandle{i+fillIdxOffset} = surf(X,Y,Z,zeros(size(X)),'FaceColor',color);
                            % remove edges and make constraint transparent
                            constraintHandle{i+fillIdxOffset}.EdgeAlpha = 0;
                            constraintHandle{i+fillIdxOffset}.FaceAlpha = 0.3;

                            case 'planar'
                            % get rectangular plane corners
                            a0 = constraints{i}.a0; % a0 ... coordinate of reference plane corner
                            a1 = constraints{i}.a1 + a0; % a1 ... first adjacent plane corner
                            a2 = constraints{i}.a2 + a0; % a2 ... second adjacent plane corner
                            a3 = a1 + a2 - a0; % a3 ... opposite plane corner
                            % generate plane coordinates
                            A = [a0 a1 a3 a2];
                            % plot plane and store fill3 handle
                            constraintHandle{i+fillIdxOffset} = fill3(A(1,:),A(2,:),A(3,:),color);
                            % remove edges and make constraint transparent
                            constraintHandle{i+fillIdxOffset}.EdgeAlpha = 0.1;
                            constraintHandle{i+fillIdxOffset}.FaceAlpha = 0.3;
                        end
                    end
                end

                % restore initial hold status
                if ~holdStatus
                    hold off
                end
                
                % draw labels
                xlabel('x'); ylabel('y'); zlabel('z');

                % disable axis limits
                set(gca,'XLim',[-Inf Inf]);
                set(gca,'YLim',[-Inf Inf]);
                set(gca,'ZLim',[-Inf Inf]);
            end
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE FUNCTIONS
    methods (Access = private)
        % setter functions
        function o = setCenterConstraint(o, fieldName, type, constraintDefinition)
            % sets constraints for coil or sensor position optimization
            %
            % INPUT:
            % fieldName - string; specifies whether the constraint is
            % active on coil- or sensor positioning with either 'coilData'
            % or 'sensorData' as input
            %
            % type - string; supported types of constraints. Currently
            % supported: 'cylindrical', 'planar'
            %
            % constraintDefinition - depends on type:
            % 'cylindrical':    (1) 3 x 1 vector; bottom center coordinate
            %                   of the cylinder
            %                   (2) 3 x 1 vector; top center coordinate of
            %                   the cylinder
            %                   (3) scalar; radius of the cylinder
            %
            %                   If no input is given, a cylinder spanning
            %                   the ROI with orientation along the z-axis
            %                   will be generated.
            %
            % 'planar':         (1) 3 x 1 vector; vector 1 spanning the
            %                   plane
            %                   (2) 3 x 1 vector; vector 2 spanning the
            %                   plane 
            %                   (3) 3 x 1 vector; plane offset from origin
            %
            %                   If no input is given, a plane slightly
            %                   below (in z-direction) the ROI with the
            %                   same dimensions in x- and y-direction will
            %                   be generated.
            
            switch type
                case 'cylindrical'
                    % check number of inputs
                    if isempty(constraintDefinition)
                        defaultShape = true;
                    elseif length(constraintDefinition) ~= 3
                        % check correct length of constraintDefinition
                        error(['Wrong number of inputs for type ''' type '''.']);
                    else
                        defaultShape = false;
                        % check correct input sizes of constraintDefinition
                        if ~all(size(constraintDefinition{1}(:))==[3 1])
                            error('Wrong size of bottom center coordinate.');
                        elseif ~all(size(constraintDefinition{2}(:))==[3 1])
                            error('Wrong size of top center coordinate.');
                        elseif length(constraintDefinition{3})~=1
                            error('Wrong size of radius.');
                        end
                    end
                                        
                    % check existence of other constraints
                    if isfield(o.(fieldName),'constraints')
                        i = length(o.(fieldName).constraints)+1;
                    else
                        i = 1;
                    end
                    
                    % if no input was given to varargin, create default
                    % cylinder
                    if defaultShape
                        % bottom center of ROI
                        constraintDefinition{1} = [mean(o.ROIData.ROI(1:2,:),2)' min(o.ROIData.ROI(3,:))] + o.ROIData.offset';
                        % top center of ROI
                        constraintDefinition{2} = [mean(o.ROIData.ROI(1:2,:),2)' max(o.ROIData.ROI(3,:))] + o.ROIData.offset';
                        % radius to ROI corners
                        constraintDefinition{3} = sqrt(sum( (o.ROIData.ROI(1:2,1) - o.ROIData.ROI(1:2,2)).^2 ))/2;
                    end
                    
                    % fill constraint fields
                    o.(fieldName).constraints{i}.type = type;
                    o.(fieldName).constraints{i}.bottomCoord = constraintDefinition{1}(:);
                    o.(fieldName).constraints{i}.topCoord = constraintDefinition{2}(:);
                    o.(fieldName).constraints{i}.radius = constraintDefinition{3};
                    
                
                case 'planar'
                    % check number of inputs
                    if isempty(constraintDefinition)
                        defaultShape = true;
                    elseif length(constraintDefinition) ~= 3
                        % check correct length of constraintDefinition
                        error(['Wrong number of inputs for type ''' type '''.']);
                    else
                        defaultShape = false;
                        % check correct input sizes of constraintDefinition
                        if ~all(size(constraintDefinition{1}(:))==[3 1])
                            error('Wrong size of first plane-spanning vector.');
                        elseif ~all(size(constraintDefinition{2}(:))==[3 1])
                            error('Wrong size of second plane-spanning vector.');
                        elseif ~all(size(constraintDefinition{3}(:))==[3 1])
                            error('Wrong size of plane offset vector.');
                        end
                    end
                                        
                    % check existence of other constraints
                    if isfield(o.(fieldName),'constraints')
                        i = length(o.(fieldName).constraints)+1;
                    else
                        i = 1;
                    end
                    
                    % if no input was given to varargin, create default
                    % plane
                    if defaultShape
                        roi = o.ROIData.ROI;
                        % first plane-spanning vector (x-length of ROI)
                        constraintDefinition{1} = [ max(roi(1,:))-min(roi(1,:));
                                        0;
                                        0];
                        % second plane-spanning vector (y-length of ROI)
                        constraintDefinition{2} = [ 0;
                                        max(roi(2,:))-min(roi(2,:));
                                        0];
                        % plane offset vector to a fifth of the ROI
                        % z-distance below the ROI
                        constraintDefinition{3} = [ min(roi(1,:));
                                        min(roi(2,:));
                                        min(roi(3,:)) - (max(roi(3,:)) - min(roi(3,:)))/5] + o.ROIData.offset(:);
                    end
                    
                    % fill constraints field
                    o.(fieldName).constraints{i}.type = type;
                    o.(fieldName).constraints{i}.a1 = constraintDefinition{1}(:);
                    o.(fieldName).constraints{i}.a2 = constraintDefinition{2}(:);
                    o.(fieldName).constraints{i}.a0 = constraintDefinition{3}(:);
                    
                    % useful additional information
                        % plane-spanning unit vectors
                    o.(fieldName).constraints{i}.a1n =  o.(fieldName).constraints{i}.a1./...
                                                    norm(o.(fieldName).constraints{i}.a1);
                    o.(fieldName).constraints{i}.a2n =  o.(fieldName).constraints{i}.a2./...
                                                    norm(o.(fieldName).constraints{i}.a2);
                    
                    
                otherwise
                    error('Type not defined.');
            end
        end
        
        % getter functions
        function o = getResolutionFromVoxels(o, voxels)
            % estimates the resolution from a given voxel grid
            % ATTENTION: works only well for voxel grids parallel to x-, y-
            % and z-directions
            %
            % INPUT:
            % voxels - 3 x nVoxels matrix; centerpoints of the voxels
            
            %%% exceptions
            % check if matrices have right sizes
            if ~any(size(voxels) == 3)
                error('Wrong size of "voxels" matrix.');
            end
            
            % transpose matrices if necessary
            if size(voxels,1) ~= 3 && size(voxels,2) == 3
                voxels = voxels';
            end
            
            %%% code
            % get resolution using the number of unique voxel coordinates
            % in x, y and z direction
            getres = @(x) length(unique(x));
            o.ROIData.resolution = [getres(voxels(1,:)) getres(voxels(2,:)) getres(voxels(3,:))];            
        end
        
%         % auxiliary functions
%         function [ newpoints, rotMatrix ] = rotatepoints(~, points, direction, startdirection)
%             % LEGACY FUNCTION -> NEWER VERSION AVAILABLE AS AUXILIARY FUNCTION
%             % rotates a pointcloud "points" with initial orientation 
%             % along the z axis (if not specified otherwise by
%             % startdirection) centered at the origin to a new directional 
%             % vector "direction".
%             % INPUT
%             % points - 3 x nPoints; pointcloud to be rotated around the
%             % origin and originally oriented along the z direction
%             %
%             % direction - 3 x 1; vector pointing to the new orientation of
%             % the pointcloud
%             %
%             % startdirection(optional) - 3 x 1; vector the rotation starts
%             % from
%             
%             %%% exceptions
%             % check input sizes
%             if ~any(size(points) == 3)
%                 error('Wrong size of "points" matrix.');
%             end
%             if ~any(size(direction) == 3)
%                 error('Wrong size of "direction" vector.');
%             end
%             
%             % transpose if necessary
%             if size(points,1) ~= 3 && size(points,2) == 3
%                 points = points';
%             end
%             direction = normcols(direction(:));
%             
%             % define start direction
%             if nargin > 3
%                 if ~any(size(startdirection) == 3)
%                     error('Wrong size of "direction" vector.');
%                 end
%                 startdirection = normcols(startdirection(:));
%             else
%                 startdirection = [0 0 1]';
%             end
%                   
%             %%% code
%             % exception: direction equals zVec
%             if all(direction==startdirection)
%                 rotMatrix = eye(3);
%                 newpoints = points;
%             elseif all(direction==-startdirection)
%                 rotMatrix = diag([1 -1 -1]);
%                 newpoints = rotMatrix*points;
%             else
%                 % get rotation axis
%                 rotAxis = normcols(cross(startdirection,direction));
%                 % get rotation angle
%                 rotAngle = acos(dot(startdirection,direction)/(norm(startdirection)*norm(direction)));
%                 % calculate rotation matrix
%                 rotMatrix = [  rotAxis(1)*rotAxis(1)*(1-cos(rotAngle))+cos(rotAngle)                    rotAxis(1)*rotAxis(2)*(1-cos(rotAngle))-rotAxis(3)*sin(rotAngle)        rotAxis(1)*rotAxis(3)*(1-cos(rotAngle))+rotAxis(2)*sin(rotAngle);
%                                 rotAxis(2)*rotAxis(1)*(1-cos(rotAngle))+rotAxis(3)*sin(rotAngle)        rotAxis(2)*rotAxis(2)*(1-cos(rotAngle))+cos(rotAngle)                    rotAxis(2)*rotAxis(3)*(1-cos(rotAngle))-rotAxis(1)*sin(rotAngle);
%                                 rotAxis(3)*rotAxis(1)*(1-cos(rotAngle))-rotAxis(2)*sin(rotAngle)        rotAxis(3)*rotAxis(2)*(1-cos(rotAngle))+rotAxis(1)*sin(rotAngle)        rotAxis(3)*rotAxis(3)*(1-cos(rotAngle))+cos(rotAngle) ];
%                 % calculate rotated points
%                 newpoints = rotMatrix*points; 
%             end
%         end
    end
end

