% clc; clear;

% initialize the toolbox
run('../../initFunctions.m')

%%%%%%%% USER INPUT:
RegionOfInterest = [0 2; 0 2; -0.1 0];
Resolution = [20 20 1];
NumberOfCoils = 20;
CoilRadius = 0.1;
NumberOfCoilSegments = 50;
NumberOfSensors = 100;
NumberOfActivationSequences = NumberOfCoils;
ActivationToBeVisualized = 1;
%%%%%%%%

%% SETUP CONFIGURATION
% initialize setup object
    % define voxel grid properties
ROI = RegionOfInterest;
res = Resolution;

    % define coil grid properties
coilCenters = rand(3,NumberOfCoils)*2;
coilCenters(3,:) = coilCenters(3,:) + CoilRadius + max(RegionOfInterest(3,:));
coilOrientations = randn(3,NumberOfCoils);
coilPattern = makeCircularCoilPattern(CoilRadius, NumberOfCoilSegments);

    % define sensor grid properties
sensorCenters = rand(3,NumberOfSensors)*2 + max(RegionOfInterest(3,:));
sensorOrientations = randn(3,NumberOfSensors);

    % configure MRXI setup
mrxisetup = SETUP();
mrxisetup = mrxisetup.setVoxelGrid(ROI, res);
mrxisetup = mrxisetup.setCoils(coilCenters, coilOrientations, coilPattern);
mrxisetup = mrxisetup.setSensors(sensorCenters, sensorOrientations);

%% ACTIVATION SEQUENCE CONFIGURATION
% calculate leadfield matrix
if NumberOfActivationSequences == NumberOfCoils
    I = eye(NumberOfCoils); % consecutive activation of each coil
else
    I = randn(NumberOfCoils,NumberOfActivationSequences); % random activations
end
[L, H_mag] = createSystemMatrix( mrxisetup, I );

% visualize setup and magnetisation
figure(1); clf;
set(gcf,'Name','MRXI setup');
mrxisetup.visualize();
hold on;
% magHandle = visualizeMagnetization( mrxisetup, H_mag, ActivationToBeVisualized );
% magHandle.Color = 'k';
daspect([1 1 1]);
hold off;

%% PHANTOM CONFIGURATION
% generate phantom data
ph1 = loadPhantom(mrxisetup.getResolution, 'Tumor');
ph2 = loadPhantom(mrxisetup.getResolution, 'P_1', 100); % SETUP.getResolution is just more convenient than SETUP.ROIData.resolution :-)
ph3 = loadPhantom(mrxisetup.getResolution, 'P_4', 100);

figure(2); clf;
set(gcf,'Name','Ground truth');
ph = {ph1, ph2, ph3};
visualizePhantom(ph);

%% CALCULATE MEASUREMENTS
% make measurements
b1 = L*ph1(:);
b2 = L*ph2(:);
b3 = L*ph3(:);
% corrupt with noise (SNR=20dB)
SNR = 20;
b1 = addGaussianNoise(b1, SNR);
b2 = addGaussianNoise(b2, SNR);
b3 = addGaussianNoise(b3, SNR);

%% RECONSTRUCTION
% normalize leadfield matrix and measurements for reconstruction
Lnorm = norm(L);
Ln = L./Lnorm;
b1n = b1./Lnorm;
b2n = b2./Lnorm;
b3n = b3./Lnorm;

% calculate inverse solution using Tikhonov regularization
alpha = 1e-5; % regularization parameter
x1 = recon_tikh(Ln, b1n, alpha);
x2 = recon_tikh(Ln, b2n, alpha);
x3 = recon_tikh(Ln, b3n, alpha);

% visualize Tikhonov reconstructions
figure(3); clf;
set(gcf,'Name','Tikhonov reconstruction');
x = {x1, x2, x3};
visualizePhantom(x, mrxisetup.getResolution);

%% calculate inverse solution using SWISTA algorithm (not a good example, but it works...)
sw1 = RECON_ISTA( Ln, b1n);
sw2 = RECON_ISTA( Ln, b2n);
sw3 = RECON_ISTA( Ln, b3n);
sw1.req.lambda = 1e-9;
sw1.req.nIter = 1.5e3;
sw1.req.nonNeg = 1;
sw1.req.deltaMin = 1e-2;
sw2.copySettings(sw1);
sw3.copySettings(sw1);

sw1.aux.xOrig = ph1;
sw1.aux.xInit = x1;
sw2.aux.xOrig = ph2;
sw2.aux.xInit = x2;
sw3.aux.xOrig = ph3;
sw3.aux.xInit = x3;

sw1.start();
sw2.start();
sw3.start();

% visualize Tikhonov reconstructions
figure(4); clf;
set(gcf,'Name','SWISTA reconstruction');
sw = {sw1.getResult(), sw2.getResult(), sw3.getResult()};
visualizePhantom(sw, mrxisetup.ROIData.resolution);
