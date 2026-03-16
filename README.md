# **mlMRXI**

## Outline

**mlMRXI** is a framework for simulating magnetorelaxometry imaging (MRXI) experiments written in MATLAB.

The framework can be used to easily configure an experimental setup including the (discretized) volume to scan, the positions, orientations, and shapes of the electromagnetic excitation coils, and the positions and orientations of the magnetometers. Multiple excitation sequences with arbitrary excitation coil currents are adjustable by the user to calculate the lead field matrix of the MRXI system. **mlMRXI** includes several predefined magnetic nanoparticle phantoms that can be used to generate artificial relaxation signals at the sensor sites.

The framework includes several reconstruction algorithms, which can be applied to recover the original magnetic particle distribution from the relaxation signals. There are also multiple convenience functions to visualize the experimental setup and the reconstruction results.

## Usage

Call

```
run('<mlMRXI_ROOT>/initFunctions.m')
```

before using **mlMRXI** in a new project. Afterwards, an MRXI experiment is simulated as follows:

- Configure your MRXI setup.

```
mrxisetup = SETUP();
mrxisetup = mrxisetup.setVoxelGrid(regionOfInterest, resolution);
mrxisetup = mrxisetup.setCoils(coilCenters, coilOrientations, coilPattern);
mrxisetup = mrxisetup.setSensors(sensorCenters, sensorOrientations);
mrxisetup.visualize();
```

- Configure your coil current activation sequences. For instance, random coil currents per activation sequence.

```
activations = randn(numberOfCoils,numberOfActivationSequences);
```

- Calculate the lead field matrix.

```
leadFieldMatrix = createSystemMatrix( mrxisetup, activations );
```

- Load a magnetic nanoparticle phantom or specify your own.

```
phantom = loadPhantom(mrxisetup.getResolution, 'Tumor');
visualizePhantom(phantom, mrxisetup.getResolution);
```

- Calculate the relaxation amplitudes measured by the magnetometers. Optionally, corrupt them with noise of your choosing.

```
relaxationAmplitudes = leadFieldMatrix * phantom(:);
% corrupt with noise (SNR=20dB)
SNR = 20;
relaxationAmplitudesNoisy = addGaussianNoise(relaxationAmplitudes, SNR);
```

- Approximate the magnetic particle phantom using one of several reconstruction algorithms from the `reconstruction/` folder. Each reconstruction class contains mandatory parameters `req` and optional parameters `aux`. All `req` parameters are initialized with suitable default values.

```
recon_obj = RECON_TIKH_ITER(leadFieldMatrix, phantom);
% adapt mandatory parameters
recon_obj.req.lambda = 1e-9;
recon_obj.req.nonNeg = 1;
% add original phantom as optional parameter to print similarity metric during reconstruction
recon_obj.aux.xOrig = phantom;
recon_obj.start();
visualizePhantom({phantom, recon_obj.getResult()}, mrxisetup.ROIData.resolution);
```

## Examples

Examples on how to construct an MRXI setup and phantoms, simulation measurements and conduct reconstructions are located under `projects/example_projects`.
