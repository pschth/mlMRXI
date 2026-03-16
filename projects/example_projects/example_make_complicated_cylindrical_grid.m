clear; clc

run('../../initFunctions.m');
%% generate PTB setup data (for sensor positions)
Resolution = [10 10 10];
mrxisetup = genPTBsetup( Resolution, 'P_3' );

%% set the voxel grid
roi = [-0.05 0.05; -0.05 0.05; -0.1 0.1];
offset = [0 0 -0.1];
orientation = [1 0 0];
mrxisetup.setVoxelGrid(roi, Resolution, offset, orientation, 'cylinder');

%% make large surrounding coil grid using different coil patterns
roidist = diff(mrxisetup.ROIData.ROI,[],2);
roicenter = mrxisetup.ROIData.ROI(:,1) + roidist./2;
ry = diff(roi(2,:))/2;
rz = diff(roi(1,:))/2;

% set parameters for coil generation
radii_y = [0.5 1 1.5 2 2.75 4.5]*0.01 + ry;
radii_z = radii_y;
coilRadii = [0.75 1 1.5 2.4 4.8 8]*0.01;
phiInitOffset = [0 pi/18 0 pi/8 pi/10 0];
radialDiscretization = [20 18 12 8 5 4];
axialDiscretization = [12 9 6 4 2 1];

% generate coils
centers = [];
orientations = [];
coilpattern = {};
coilpatternassignment = [];
for p = 1:length(radii_y)
    phi = linspace(0,2*pi,radialDiscretization(p)+1)+phiInitOffset(p);
    phi = phi(1:end-1);
    x = linspace(mrxisetup.ROIData.ROI(3,1), mrxisetup.ROIData.ROI(3,2), axialDiscretization(p)*2+1);
    x = x(2:2:end-1);
    y = radii_y(p) * cos(phi);
    z = radii_z(p) * sin(phi);
    
    if p == length(radii_y)
        y(2) = [];
        z(2) = [];
    end

    X = kron(x, ones(1, length(y)));
    Y = repmat(y, 1, length(x));
    Z = repmat(z, 1, length(x));
    orientations = [orientations [zeros(size(X(:))), -Y(:), -Z(:)]']; %#ok<*AGROW>
    Z = Z + offset(3);
    centers = [centers [X; Y; Z]];
    coilpattern = [coilpattern {makeSpiralCoilPattern(coilRadii(p), [], 5e-4)}];
    coilpatternassignment = [coilpatternassignment p*ones(1,length(X))];
end

% add coils to MRXI setup object
mrxisetup.setCoils(centers, orientations, coilpattern, coilpatternassignment);


%% visualize the setup
f1 = figure(1);clf;
f1.Color = ones(1,3);
mrxisetup.visualize;
daspect([1 1 1]);
view([-81 11]);
drawnow;
