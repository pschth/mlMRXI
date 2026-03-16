function [voxHandle, coilHandle, sensHandle] = visualizeSetup( SETUP, VoxAsScatter )
% visualizes current setup configuration using gray voxels, red coils, and
% blue sensors
% 
% INPUT:
% SETUP - MRXI setup configuration object
% 
% VoxAsScatter (optional) - boolean; if true, voxels are represented as
% scatter plot rather than the ROI as box; default: false
% 
% OUTPUT:
% voxHandle - scatter or patch handle, depending on if voxels are displayed
% as points or cuboids
% 
% coilHandle - nCoils x 1 cell; contains the plot handles for all drawn
% coils
% 
% sensorHandle - 1 x 2 cell; contains the scatter and the quiver plot
% handles for the drawn sensors

    % use cuboids to visualize voxels by default
    if nargin < 2
        VoxAsScatter = false;
    elseif isempty(VoxAsScatter)
        VoxAsScatter = false;
    end
    
    % set to scatter visualization if voxel dimensions are not provided
    if ~isfield(SETUP.ROIData, 'voxelDims')
        VoxAsScatter = true;
    end
    
    % get initial hold status
    holdflag = ishold;
    % clear figure if holdflag is FALSE
    if ~holdflag
        clf;
    end
    
    if SETUP.getNumberOfVoxels == 0
        % if no voxels are defined, return empty voxHandle
        voxHandle = [];
    else
        % get actual ROI position
% % % % % to catch exceptions from old SETUP objects. INFO FOR pyMRXI: not
% % % % % necessary!
        if ~isfield(SETUP.ROIData, 'rotvector')
            SETUP.ROIData.rotvector = [0 0 1]';
        end
        if ~isfield(SETUP.ROIData, 'offset')
            SETUP.ROIData.offset = [0 0 0]';
        end
% % % % %         
                    
        if VoxAsScatter
            % represent voxels as scatter plot if demanded
            voxHandle = scatter3(SETUP.ROIData.voxels(1,:), ...
                                 SETUP.ROIData.voxels(2,:), ...
                                 SETUP.ROIData.voxels(3,:), ...
                                 30, 'filled', 'k');
            % make gray scatter points with black edges
            voxHandle.MarkerEdgeColor = [0 0 0];
            voxHandle.MarkerFaceColor = [0.75 0.75 0.75];
        elseif ~isfield(SETUP.ROIData, 'shape') || ...
                 strcmp(SETUP.ROIData.shape, 'cuboid') || ...
                 strcmp(SETUP.ROIData.shape, 'cylinder') || ...
                 strcmp(SETUP.ROIData.shape, 'sphere')
            % represent voxel hull as patches
            voxHandle = drawvoxels(SETUP);
        end
        hold on
    end

%     visualize coils if available
    if SETUP.getNumberOfCoils == 0
        coilHandle = [];
    else
        coilHandle = visualizeCoils(SETUP);
        hold on;
    end
    
%     visualize sensors if available
    if SETUP.getNumberOfSensors == 0
        sensHandle = [];
    else
        sensHandle = visualizeSensors( SETUP, 'Color', 'b');
    end
    
%     make legend entries for nonempty voxel, coil, and sensor handles
    legendHandles = [];
    legendEntries = {};
    % add voxel legend, if necessary
    if ~isempty(voxHandle)
        legendHandles = voxHandle;
        legendEntries = {'Voxels'};
    end
    % add coil legend, if necessary
    if ~isempty(coilHandle)
        legendHandles = [legendHandles coilHandle{ceil(end/2)}];
        legendEntries = [legendEntries {'Coils'}];
    end
    % add sensor legend, if necessary
    if ~isempty(sensHandle)
        legendHandles = [legendHandles sensHandle{1}];
        legendEntries = [legendEntries {'Sensors'}];
    end
    % display legend
    if ~isempty(legendHandles)
        legend(legendHandles, legendEntries);
    end
    
    % reset initial hold status
    if ~holdflag
        hold off
    end

    % set tight limits to increase visibility
    axis tight
    % turn view of setup to what is usually not bad...
    view([40 25])
end

