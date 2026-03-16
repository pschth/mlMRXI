function figHandle = visualizeMagnetization( mrxisetup, H, activation, normalizeArrows, arrowScale )
% visualizes the magnetic field vectors 'H' in the voxels given in
% 'mrxisetup'
% 
% INPUT:
% mrxisetup - MRXI setup configuration object
% 
% H - 3 x nVoxels*nActivations matrix; contains magnetization vectors of
% each voxel during each activation
% 
% activation - scalar (optional); number of the activation sequence to be
% visualized; default: 1
% 
% normalizeArrows - boolean (optional); if true, magnetization vector
% lengths are normalized to unit length; default: false
% 
% arrowScale - scalar (optional); arrow scaling factor; default: 1
% 
% OUTPUT:
% figHandle - plot handle of the quiver plot

    if nargin < 3
        activation = 1;
    elseif isempty(activation)
        activation = 1;
    end
    
    if nargin < 4
        normalizeArrows = false;
    elseif isempty(normalizeArrows)
        normalizeArrows = false;
    end
    
    if nargin < 5
        arrowScale = 1;
    elseif isempty(arrowScale)
        arrowScale = 1;
    end
        
    
%     transpose H if necessary
    if size(H,1) ~= 3 && size(H,2) == 3
        H = H';
    end
    
%     extract desired activation sequence from H
    nVox = size(mrxisetup.ROIData.voxels,2);
    H = H(:, (activation-1)*nVox+1:activation*nVox);
    
%     normalize vectors if demanded
    if normalizeArrows
        H = normcols(H);
    end
    
%     visualize magnetization
    figHandle = quiver3(mrxisetup.ROIData.voxels(1,:),mrxisetup.ROIData.voxels(2,:),mrxisetup.ROIData.voxels(3,:),...
        H(1,:),H(2,:),H(3,:), arrowScale);
end

