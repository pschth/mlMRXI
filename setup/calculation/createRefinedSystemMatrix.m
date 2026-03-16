function [L_refined, L_refined_whole] = ...
    createRefinedSystemMatrix(refinementLvl, SETUP, I, thresh_RAM, verbose )
% calculates a refined version "L_refined" of the MRX forward model system
% matrix. Only works for SETUPs with ROIData defined via the function
% "setVoxelGrid".
% 
% refinementLvl: integer, defines the level of refinement. E.g. 2 means a
% single voxel in the original SETUP is represented in L_refined as mean
% value of a 2x2x2 grid instead of only the voxel center
% 
% SETUP: standard MRXI config file from the MRXI toolbox
% 
% I: nCoils x nActivations matrix, contains coil current values per
% activation. Default: Identity-matrix
% 
% thresh_RAM: a RAM-threshold in Gb can be stated. If the
% predicted RAM-usage of the calculation is below, the fast version of the
% algorithm is computed. Otherwise, the slower, RAM-efficient version is
% computed. Default: 4 Gb
% 
% verbose: boolean; if false, no text messages are generated. Default: true
% 

    if nargin < 3
        I = eye(SETUP.getNumberOfCoils);
    elseif isempty(I)
        I = eye(SETUP.getNumberOfCoils);
    end

    if nargin < 4
        thresh_RAM = 4;
    elseif isempty(thresh_RAM)
        thresh_RAM = 4;
    end
    
    if nargin < 5
        verbose=true;
    elseif isempty(verbose)
        verbose=true;
    end
    
    % check if refinementLvl is integer
    if refinementLvl < 1 || refinementLvl/floor(refinementLvl) ~= 1
        error('refinementLvl is not integer or smaller than 1.');
    end
    
    if refinementLvl == 1
        % if refinementLvl = 1 calculate normal system matrix and leave
        % refined system matrix blank
        L_refined_whole = [];
        L_refined = createSystemMatrix(SETUP, I, thresh_RAM, verbose);
    else
        % copy setup
        mrxi = SETUP.copySetup;
        nVoxelOrig = mrxi.getNumberOfVoxels;
        % get voxel dimensions
        resOrig = mrxi.getResolution;
        voxDim = diff(mrxi.ROIData.ROI,[],2)./resOrig';
        % initialize new voxel matrix
        multRefLvl = refinementLvl^3;
        voxelsRef = zeros(3, nVoxelOrig * multRefLvl);
        % compute refined voxel positions
        for vOrig = 1:mrxi.getNumberOfVoxels
            % get voxel borders
            voxMin = mrxi.ROIData.voxels(:,vOrig) - voxDim./2;
            voxMax = mrxi.ROIData.voxels(:,vOrig) + voxDim./2;
            % get refined voxel coordinates
            xRef = linspace(voxMin(1), voxMax(1), refinementLvl + 2); xRef = xRef(2:end-1);
            yRef = linspace(voxMin(2), voxMax(2), refinementLvl + 2); yRef = yRef(2:end-1);
            zRef = linspace(voxMin(3), voxMax(3), refinementLvl + 2); zRef = zRef(2:end-1);
            % mesh grid with new coordinates
            [X,Y,Z] = meshgrid(xRef, yRef, zRef);
            % store refined voxel positions
            vRef = (vOrig-1)*multRefLvl+1 : vOrig*multRefLvl;
            voxelsRef(:,vRef) = [X(:) Y(:) Z(:)]';
        end

        % set refined voxel coordinates
        mrxi.setVoxelCoordinates(voxelsRef);
        % get refined leadfield matrix
        L_refined_whole = createSystemMatrix(mrxi, I, thresh_RAM, verbose );
        % average over refined voxels
        L_refined = squeeze(mean( reshape(L_refined_whole, ...
            size(L_refined_whole,1), multRefLvl, nVoxelOrig), 2 ));
    end
end
