function out = getVoxelSensorRelation( SETUP, verbose )
% calculates the vectorial dipole relation of each voxel to each sensor
% 
% INPUT:
% SETUP - MRXI setup configuration object
% 
% verbose (optional) - boolean; if false, no text messages are generated. Default: true
% 
%
% OUTPUT:
% out - 3 x nSensors x nVoxels matrix; vectorial dipole relation of each
% voxel to each sensor
%    
    if nargin < 2
        verbose = true;
    elseif isempty(verbose)
        verbose = true;
    end
        
    % get calculation constants
    nVoxels = SETUP.getNumberOfVoxels;
    nSensors = SETUP.getNumberOfSensors;

    % initialize all huge matrices to predict RAM usage
    out = zeros([3 nVoxels nSensors]);
         
    if verbose
        fprintf('\n\tCalculating theoretical influences of each sensor on each voxel... ');
    end
    % for all coils
    for si = 1:nSensors
        % get all voxel-coilCenter distance vectors
        vsd = bsxfun(@minus, SETUP.sensorData.centers(:,si), SETUP.ROIData.voxels);

        % get voxel-coilCenter distances
        vsdNorm = sqrt(sum(vsd.*vsd,1));

        % for each coil get theoretical influences of each coil on each voxel
        out(:,:,si) = bsxfun(@times,3*vsd,(SETUP.sensorData.orientations(:,si)'*vsd)./vsdNorm.^5)-...
                       kron(SETUP.sensorData.orientations(:,si), 1./vsdNorm.^3);       
    end

    mu_0 = 4*pi * 1e-7;
    out = permute(out, [1 3 2])*mu_0./(4*pi);
    
    if verbose
        fprintf('\nDone!\n');
    end
end