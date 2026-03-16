function [phantom, voxelDelIdx] = loadPhantom( Resolution, phantomType, MNPQuantity, shape )
% creates/loads predefined MNP phantoms
% 
% INPUT:
% Resolution - 3 x 1 vector; voxel grid resolution
% 
% phantomType - string; defines the phantom type to load/create. Valid
% strings: 
% 'P_1' to 'P_5' ... PTB's P-phantom from bottom z-plane (P_1) to top
% z-plane (P_5)
% 'PTB' or 'PTB_boost' ... PTB's PTB-phantom with 'P' in the top z-plane,
% 'T' in the central z-plane and 'B' in the bottom z-plane
% 'SL' ... Shepp-Logan phantom
% 'MSL' ... Modified Shepp-Logan phantom with improved contrast
% 'Tumor' ... Tumor phantom defined in collaboration with Janic Föcke,
% University of Münster
% 'Checkerboard' ... Checkerboard phantom with 5x5x5 grid
% 'Planes' ... Planes phantom with alternating planes in 5 z planes. Can be
% extended with '_x', '_y' or '_z' (default) to rotate towards the
% respective direction. 
% 
% MNPQuantity (optional) - specifies the overall amount of MNPs inside the
% phantom in [mg]; default: 1 [mg]; ATTENTION: output 'phantom' in [kg]
% 
% shape (optional) - can be 'cuboid' or 'cylinder'. Default: 'cuboid'
% 
% OUTPUT:
% phantom - Resolution sized matrix, contains the phantom data in [kg]
% 
% voxelDelIdx - numel(Resolution) sized boolean vector, contains the voxels
% marked for deletion if using a cylindrical phantom shape to go from the
% Resolution sized, cuboidal grid to the cylindrical grid. Required for
% phantom visualization.

%%% exceptions
if nargin < 4 || isempty(shape)
    shape = 'cuboid';
end
if ~(strcmp(shape,'cuboid') || strcmp(shape,'cylinder') || strcmp(shape,'sphere'))
    error('Shape not defined.');
end
            
if nargin < 3 || isempty(MNPQuantity)
    MNPQuantity = 1;
end

if ~all(size(Resolution(:))==[3 1])
    error('No valid input type of Resolution. Use 3 x 1 integer vector.');
end

%%% code

phantomTypeSplit = strsplit(phantomType,'_');

% create initial phantom data
switch phantomTypeSplit{1}
    case 'P'
        % PTB P-phantom
        phantom = zeros(100,100,100);
        
        if length(phantomTypeSplit)~=2
            error('P-phantom has to be P_1, P_2, P_3, P_4, P_5, or P_32.');
        end
        PIdx = str2double(phantomTypeSplit{2});
        if  ~(PIdx>=1 && PIdx<=5 || PIdx==32)
            error('P-phantom has to be P_1, P_2, P_3, P_4, P_5, or P_32.');
        end
        
        if PIdx>=1 && PIdx<=5
            z = PIdx*20-19:PIdx*20;
            phantom(61:70,21:90,z) = 1;
            phantom(31:40,21:60,z) = 1;
            phantom(41:60,21:30,z) = 1;
            phantom(41:60,51:60,z) = 1;
        else
            z = 3*20-19:3*20;
            phantom(61:70,21:80,z) = 1;
            phantom(31:40,21:60,z) = 1;
            phantom(41:60,21:30,z) = 1;
            phantom(41:60,51:60,z) = 1;
            phantom = permute(phantom, [2 1 3]);
            phantom = phantom(:,end:-1:1,:);
        end
        
    case 'PTB'
        % PTB P-phantom
        phantom = zeros(100,100,100);
        
        if all(length(phantomTypeSplit)~=[1 2])
            error('PTB-phantom has to be either PTB or PTB_boost.');
        end
        if length(phantomTypeSplit) == 2
            if  ~strcmp(phantomTypeSplit{2}, 'boost')
                error('PTB-phantom has to be either PTB or PTB_boost.');
            end
        end
        % P
        z = 5*20-19:5*20;
        phantom(11:20,11:80,z) = 1;
        phantom(21:50,41:50,z) = 1;
        phantom(21:50,71:80,z) = 1;
        phantom(41:50,51:70,z) = 1;
        % T
        z = 3*20-19:3*20;
        phantom(41:50,11:80,z) = 1;
        phantom(11:80,81:90,z) = 1;
        % B
        z = 1*20-19:1*20;
        phantom(51:90,11:20,z) = 1;
        phantom(51:90,41:50,z) = 1;
        phantom(51:90,71:80,z) = 1;
        phantom(51:60,11:80,z) = 1;
        phantom(81:90,11:80,z) = 1;
        % flip to correct orientation
        phantom = phantom(end:-1:1, end:-1:1, :);
       
    case 'SL'
        % Shepp-Logan phantom
        phantom = phantom3d('Shepp-Logan',max([Resolution 100]));
        
    case 'MSL'
        % Modified Shepp-Logan phantom with improved contrast
        phantom = phantom3d('Modified Shepp-Logan',max([Resolution 100]));
        
    case 'Tumor'
        % Tumor phantom
        e =    [...
				1  .6900  .620  .900      0       0       0      0      0      0; ...
				-1  1	0.1  .2      -.7       0       0      30      0      0; ...
				-1  1	0.1  .2      1.0       -.1       -.1      0       30     50; ...
				-1  1	0.1  .2      1.0       -.1       -.1      0       30     -50; ...
				];
			
        phantom = phantom3d(e, max([Resolution 100]));
        phantom(phantom < 0) = 0;
        
    case 'Checkerboard'
        % Checkerboard phantom with 5x5x5 grid
        phantom = zeros(100,100,100);
        
        z1 = [1:20 41:60 81:100];
        z2 = setdiff(1:100, z1);
        phantom(z1, z1, z1) = 1;
        phantom(z2, z2, z2) = 1;
        
    case 'Planes'
        
        % Planes phantom with alternating planes in 5 z planes
        phantom = zeros(100,100,100);

        z1 = [1:20 41:60 81:100];
        z2 = setdiff(1:100, z1);
        x = 21:80;
        y1 = 11:50;
        y2 = 51:90;
        phantom(x, y1, z1) = 1;
        phantom(x, y2, z2) = 1;    
        
        if length(phantomTypeSplit) > 1
            switch phantomTypeSplit{2}
                case 'x'
                    phantom = permute(phantom, [3 2 1]);
                case 'y'
                    phantom = permute(phantom, [2 3 1]);
                case 'z'
                otherwise
                    error('Planes phantom can be only be oriented along ''x'', ''y'' and ''z''.');
            end
        end 
        
    case 'ConcentrationSeries'
        
        % Concentration series phantom from PTB measurements
        
        % first create 2D high-res plane of MNP distribution
        Res2D = 1001; % 2D resolution
        phantomFlat = zeros(Res2D); % xy-plane of ROI (12cm x 12cm)
        l = 0.12; % square phantom container dimension
        r = 0.004; % inner container radius
        rd = r/l*Res2D; % discretized radius
        MNPamount = [2.5 2 0 1 0.5 0.25]; % mg MNP per container
        
        [X,Y] = meshgrid(1:Res2D, 1:Res2D);
        
        step = Res2D/7; % distance between MNP containers
        
        i = 1;
        for y = 6*step:-step:step
            c = [round(Res2D/2) y]; % center of container
            distToC = reshape(...
                        sqrt(sum(([X(:) Y(:)] - repmat(c,Res2D^2,1)).^2,2)),...
                        Res2D, Res2D);
            phantomFlat(distToC<=rd) = MNPamount(i);
            i=i+1;
        end
        
        % interpolate resolution to smaller grid and extrude
        Res3D = 101;
        [X,Y] = meshgrid(linspace(1,Res2D,Res3D), linspace(1,Res2D,Res3D));
        phantomFlat = interp2(phantomFlat, X,Y,'linear');
        
        phantom = zeros(Res3D, Res3D, Res3D);
        zd = [round(0.026/0.06*Res3D) round(0.036/0.06*Res3D)]; % discretized z coordinates of phantom containers
        phantom(:,:,zd(1):zd(2)) = repmat(phantomFlat,1,1,zd(2)-zd(1)+1);
        
    case 'VolumeSeries'
        
        % Volume series phantom from PTB measurements
        
        % first create 2D high-res plane of MNP distribution
        Res2D = 1001; % 2D resolution
        phantomFlat = zeros(Res2D); % xy-plane of ROI (12cm x 12cm)
        l = 0.12; % square phantom container dimension
        r = [5 4 3 2 1.5 1] * 1e-3; % inner container radii
        rd = r./l*Res2D; % discretized radius
        
        [X,Y] = meshgrid(1:Res2D, 1:Res2D);
        
        step = Res2D/7; % distance between MNP containers
        
        i = 1;
        for y = 6*step:-step:step
            c = [round(Res2D/2) y]; % center of container
            distToC = reshape(...
                        sqrt(sum(([X(:) Y(:)] - repmat(c,Res2D^2,1)).^2,2)),...
                        Res2D, Res2D);
            phantomFlat(distToC<=rd(i)) = 1;
            i=i+1;
        end
        
        % interpolate resolution to smaller grid and extrude
        Res3D = 101;
        [X,Y] = meshgrid(linspace(1,Res2D,Res3D), linspace(1,Res2D,Res3D));
        phantomFlat = interp2(phantomFlat, X,Y,'linear');
        
        phantom = zeros(Res3D, Res3D, Res3D);
        zd = [round(0.026/0.06*Res3D) round(0.036/0.06*Res3D)]; % discretized z coordinates of phantom containers
        phantom(:,:,zd(1):zd(2)) = repmat(phantomFlat,1,1,zd(2)-zd(1)+1);
        
        
    case 'Cross'
        
        % central MNP columns in every direction
        
        phantom = zeros(99,99,99);
        phantom(45:55, 45:55, :) = 1;
        phantom(45:55, :, 45:55) = 1;
        phantom(:, 45:55, 45:55) = 1;
        
    case 'Sphere'
        % Sphere-phantom
        % Usage: 'Sphere_xCoord_yCoord_zCoord_radius'
        phantom = zeros(100,100,100);
        
        if length(phantomTypeSplit)~=2
            error('P-phantom has to be P_1, P_2, P_3, P_4, P_5, or P_32.');
        end
        PIdx = str2double(phantomTypeSplit{2});
        if  ~(PIdx>=1 && PIdx<=5 || PIdx==32)
            error('P-phantom has to be P_1, P_2, P_3, P_4, P_5, or P_32.');
        end
        
        if PIdx>=1 && PIdx<=5
            z = PIdx*20-19:PIdx*20;
            phantom(61:70,21:90,z) = 1;
            phantom(31:40,21:60,z) = 1;
            phantom(41:60,21:30,z) = 1;
            phantom(41:60,51:60,z) = 1;
        else
            z = 3*20-19:3*20;
            phantom(61:70,21:80,z) = 1;
            phantom(31:40,21:60,z) = 1;
            phantom(41:60,21:30,z) = 1;
            phantom(41:60,51:60,z) = 1;
            phantom = permute(phantom, [2 1 3]);
            phantom = phantom(:,end:-1:1,:);
        end
        
    otherwise
        
        error(['PhantomType ' phantomType ' not supported.']);
        
end

% interpolate initial phantom to actual grid size
phantom = interpolatePhantomToNewResolution(phantom, size(phantom), Resolution);

% cut out cylinder if specified shape is 'cylinder'
% calculate voxel coordinates
% first, calculate voxels with relative positions identical to the
% cylindrical voxel placement in SETUP with radius=1
x = linspace(-1,1,Resolution(1)*2+1);
x = x(2:2:end-1);
y = linspace(-1,1,Resolution(2)*2+1);
y = y(2:2:end-1);
% make grid from coordinates
if strcmp(shape,'cylinder')
    % if shape is cylinder, remove all voxels outside of cylinder
    % radius along z-axis
    [X,Y] = meshgrid(y,x);
    R = sqrt(X.^2 + Y.^2); % calculate radius to central axis at every point
    voxelDelIdx = repmat(R>1, 1, 1, Resolution(3)); % mark every coordinate greater 1 for deletion
    phantom = phantom(~voxelDelIdx);    
else
    voxelDelIdx = [];
end

% rescale to absolute MNP amount inside phantom
phantom = phantom./(sum(phantom(:))*1e6).*MNPQuantity;
phantom(isnan(phantom)) = 0;

end
