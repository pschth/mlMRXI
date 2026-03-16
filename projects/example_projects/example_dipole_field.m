clc; clear;

% initialize the toolbox
run('../../initFunctions.m')

%%%%%%%% USER INPUT:
RegionOfInterest = [0 2; 0 2; -2 0];
Resolution = [10 10 10];
NumberOfSensors = 1;
%%%%%%%%

%% SETUP CONFIGURATION
% initialize setup object
    % define voxel grid properties
ROI = RegionOfInterest;
res = Resolution;

    % define coil grid properties
coilCenters = [mean(ROI(1:2,:),2); ROI(3,2)+0.1];
coilOrientations = [0;0;-1];

    % define sensor grid properties (sensors are actually not required for
    % this example, but nonetheless defined for demonstration of the
    % SETUP() configuration)
sensorCenters = rand(3,NumberOfSensors)*2 + max(RegionOfInterest(3,:));
sensorOrientations = randn(3,NumberOfSensors);

    % configure MRXI setup
mrxisetup = SETUP();
mrxisetup = mrxisetup.setVoxelGrid(ROI, res);
mrxisetup = mrxisetup.setCoils(coilCenters, coilOrientations); % .setCoils without definition of one or multiple coil patterns (3rd argument) sets dipolar sources instead of coils
mrxisetup = mrxisetup.setSensors(sensorCenters, sensorOrientations);

%% ACTIVATION SEQUENCE CONFIGURATION
% calculate leadfield matrix
H_mag = getMagneticField(mrxisetup);

% visualize setup and magnetisation
figure(1); clf;
set(gcf,'Name','MRXI setup');
[voxHandle, coilHandle, sensHandle] = mrxisetup.visualize();
voxHandle.FaceAlpha = 0.2;
voxHandle.EdgeAlpha = 0.1;
hold on;
magHandle = visualizeMagnetization( mrxisetup, H_mag, 1, true ); % set fourth argument to false for correct vector lengths
magHandle.Color = 'k';
daspect([1 1 1]);
hold off;
title('Visualization of the dipole''s magnetic field (normalized vectors)')