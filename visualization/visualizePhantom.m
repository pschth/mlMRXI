function visualizePhantom( phantom, Resolution, UseCommonCRange, voxelDelIdx )
% visualizes MNP distributions 
% x-Axis (first dimension, ascending) horizontal, left to right
% y-Axis (second dimension, ascending) vertical, top to bottom
% z-Axis (first dimension, ascending) vertical images, left to right
% 
% INPUT:
% phantom - cell; contains an arbitrary number of 3D matrix (one per cell)
% OF THE SAME SIZES containing the MNP amount per voxel, respectively.
% 
% Resolution (optional) - 3 x 1 vector; voxel grid resolution
%
% UseCommonCRange (optional) - bool; defines if color range of phantoms is
% defined by each phantom individually (false, Default) or commonly (true)
% 
% voxelDelIdx - numel(Resolution) sized boolean vector, contains the voxels
% marked for deletion if using a phantom shape other than cuboidal to go
% from the Resolution sized, cuboidal grid to the desired grid.

if nargin < 4
    delVoxels = false;
elseif isempty(voxelDelIdx)
    delVoxels = false;
elseif numel(voxelDelIdx) == prod(Resolution)
    delVoxels = true;
else
    error('voxelDelIdx has not the dimensions specified in Resolution');
end
    
            
if nargin < 3
    UseCommonCRange = false;
elseif isempty(UseCommonCRange)
    UseCommonCRange = false;
end


if nargin > 1 && ~isempty(Resolution)
    if ~all(size(Resolution(:))==[3 1])
        error('Resolution wrong size.');
    end
    doReshape = true;
else
    doReshape = false;
end

if iscell(phantom)
    nPh = length(phantom);
    if delVoxels
        for ph = 1:nPh
            tempPhant = zeros(Resolution);
            tempPhant(voxelDelIdx) = nan;
            tempPhant(~voxelDelIdx) = phantom{ph};
            phantom{ph} = tempPhant;
        end
    elseif doReshape
        for ph = 1:nPh
            if numel(size(phantom{ph})) ~= numel(Resolution)
                phantom{ph} = reshape(phantom{ph}, Resolution);
            elseif ~all(size(phantom{ph}) == Resolution)
                phantom{ph} = reshape(phantom{ph}, Resolution);
            end
        end
    end
    nZ = size(phantom{1},3);
    
    % get common color range if desired
    if UseCommonCRange
        cRange = [Inf, -Inf];
        for ph = 1:nPh
            currRange = [min(phantom{ph}(:),[],'omitnan') max(phantom{ph}(:),[],'omitnan')];
            if cRange(1) > currRange(1)
                cRange(1) = currRange(1);
            end
            if cRange(2) < currRange(2)
                cRange(2) = currRange(2);
            end
        end
        if delVoxels
            delVoxValue =  cRange(1) - diff(cRange)*257.1/256^2;
            cRange(1) = delVoxValue;
        end
    end
    
    % plot z-planes of phantoms
    for ph = 1:nPh
        if ~UseCommonCRange
            cRange = [min(phantom{ph}(:),[],'omitnan') max(phantom{ph}(:),[],'omitnan')];
            if delVoxels
                delVoxValue =  cRange(1) - diff(cRange)*257.1/256^2;
                cRange(1) = delVoxValue;
            end
        end
        if delVoxels
            phantom{ph}(voxelDelIdx) = delVoxValue;
            colormap([0.91 0.91 0.91; parula(256)]);
        else
            colormap(parula(256));
        end
        for z = 1:nZ
            subaxis(nPh, nZ, z, ph, ...
                'SH', 0, 'SV', 0.005,...
                'MB', 0, 'MT', 0, 'MR', 0, 'ML', 0);
            imagesc(phantom{ph}(:,:,z)');
            legend('hide');
            set(gca, 'XTick', [], 'XTickLabel', [],...
                     'YTick', [], 'YTickLabel', []);
            if cRange(1) ==  cRange(2)
                cRange(2) = cRange(1) + eps;
            end
            caxis(cRange);
        end
    end
    
else
    % plot z-planes of phantom
    nPh = 1;
    cRange = [min(phantom(:),[],'omitnan') max(phantom(:),[],'omitnan')];
    if delVoxels
        delVoxValue =  cRange(1) - diff(cRange)*257/256^2;
        cRange(1) = delVoxValue;
        
        tempPhant = zeros(Resolution);
        tempPhant(voxelDelIdx) = delVoxValue;
        tempPhant(~voxelDelIdx) = phantom;
        phantom = tempPhant;
        colormap([1 1 1;parula(256)]);
    elseif doReshape
        phantom = reshape(phantom, Resolution);
        colormap(parula(256));
    end
    nZ = size(phantom,3);
    for z = 1:nZ
        subaxis(1, nZ, z,...
                'SH', 0.005, 'SV', 0,...
                'MB', 0, 'MT', 0, 'MR', 0, 'ML', 0);
        imagesc(phantom(:,:,z)');
        legend('hide');
        set(gca, 'XTick', [], 'XTickLabel', [],...
                 'YTick', [], 'YTickLabel', []);
        if cRange(2) > cRange(1)
            caxis(cRange);
        end
        
    end
end

monDims = get(0,'MonitorPositions'); % get monitor dimensions
maxHeight = min(monDims(:,4))-100; % get lowest display height
maxWidth = min(monDims(:,3)); % get lowest display width

% define good/maximum dimensions of figure
axDim = 0.95 * min(maxHeight/nPh, maxWidth/nZ);
figWidth = axDim*nZ;
figHeight = axDim*nPh+0.005*axDim*(nPh-1);
figPos = get(gcf,'Position');

% set final figure dimensions
set(gcf,'Position', [figPos(1:2) figWidth figHeight]);





end

