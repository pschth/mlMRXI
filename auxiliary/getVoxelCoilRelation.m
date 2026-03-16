function out = getVoxelCoilRelation( SETUP, verbose )
% calculates the vectorial relation of each voxel to each coil independent
% of coil currents (i.e. all coil currents = 1).
% 
% INPUT:
% SETUP - MRXI setup configuration object
% 
% verbose (optional) - boolean; if false, no text messages are generated. Default: true
% 
%
% OUTPUT:
% out - 3 x nCoils x nVoxels matrix; vectorial relation of each
% voxel to each coil


    if nargin < 2
        verbose = true;
    elseif isempty(verbose)
        verbose = true;
    end
    % get magnetic field vectors of all coils in all voxels in a
    % nVoxel*nCoil x 3 matrix
    Hmag = getMagneticField(SETUP,[],verbose);
    
    % reshape matrix to 3 x nCoils x nVoxels
    nVoxel = SETUP.getNumberOfVoxels;
    nCoils = SETUP.getNumberOfCoils;
    out = permute(reshape(Hmag', 3, nVoxel, nCoils), [1 3 2]);
end