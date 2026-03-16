function figHandle = drawvoxels(SETUP)
% draws the outer hull of a voxel grid as gray patches with black lines
% connecting adjacent voxel corners
% 
% INPUT:
% SETUP - MRXI setup configuration object


    % get minimum and maximum voxel corner points
    l = -SETUP.ROIData.voxelDims./2; % minimum corner point
    u = SETUP.ROIData.voxelDims./2; % maximum corner point
    
    % define patch coordinates centered at origin: per coordinate 6 rows (6
    % faces of the cuboidal voxel) and 4 rows (4 corners per face)
    X = [l(1) l(1) l(1) l(1); l(1) l(1) u(1) u(1); l(1) l(1) u(1) u(1);...
         u(1) u(1) u(1) u(1); l(1) l(1) u(1) u(1); l(1) l(1) u(1) u(1)];
    Y = [l(2) l(2) u(2) u(2); l(2) l(2) l(2) l(2); l(2) u(2) u(2) l(2);...
         l(2) l(2) u(2) u(2); u(2) u(2) u(2) u(2); l(2) u(2) u(2) l(2)];
    Z = [l(3) u(3) u(3) l(3); l(3) u(3) u(3) l(3); l(3) l(3) l(3) l(3);...
         l(3) u(3) u(3) l(3); l(3) u(3) u(3) l(3); u(3) u(3) u(3) u(3)];
     
    % rotate patch coordinates to orientation of the ROI
    pCoords = rotatepoints([X(:) Y(:) Z(:)], SETUP.ROIData.rotvector);
    
    % get number of voxels
    nVoxels = SETUP.getNumberOfVoxels;
    % initialize patch coordinates (like above) for all voxels in a patch
    % coordinate collection
    XAll = zeros(4, 6*nVoxels);
    YAll = XAll;
    ZAll = XAll;
    % for every voxel, translate patch coordinates to actual position and
    % draw box for each voxel
    for v = 1:nVoxels
        fillIdx = (v-1)*6+1:v*6; % index for filling the matrices
        
        % translate centered patch coordinates to actual voxel position v
        pCoordsCurr = bsxfun(@plus, pCoords, SETUP.ROIData.voxels(:,v));
        
        % insert translated patch coordinates into patch coordinate
        % collection
        XAll(:, fillIdx) = reshape(pCoordsCurr(1,:), 6, 4)';
        YAll(:, fillIdx) = reshape(pCoordsCurr(2,:), 6, 4)';
        ZAll(:, fillIdx) = reshape(pCoordsCurr(3,:), 6, 4)';
    end
    % draw ALL voxel faces as patches
    figHandle = patch('XData',XAll, 'YData',YAll, 'ZData', ZAll);

% % %     remove INNER patches (all patches containing vertices with 24
% % %     connections)
    % assign identical indices to identical vertex coordinates (=voxel
    % corner points)
    [~, ~, ic] = unique( ... 
                    round(figHandle.Vertices, 9), ... round coordinates to 9th decimal place (=nanometers) to remove numerical imprecisions that might disturb 'unique'
                    'rows');
    % mark every vertex index that occurs exactly 24 times (3 patches in 8
    % quadrants for every inner corner point) TRUE for deletion (=inner
    % vertices of the voxel grid)
    delIdxUniqueVertBool = histcounts(ic,'BinMethod','integers') == 24;
    % get decimal indices of TRUE boolean indices for unique vertices to
    % delete
    delIdxUniqueVert = find(delIdxUniqueVertBool);
    % mark every vertex that is part of delIdxUniqueVert (=24 occurences)
    % for deletion
    delIdxVertBool = ismember(ic, delIdxUniqueVert);
    % get decimal indices of TRUE boolean indices for vertices to delete
    delIdxVert = find(delIdxVertBool);
    % delete all faces (=rows) that contain at least one of the marked
    % vertices
    figHandle.Faces( any(ismember(figHandle.Faces, delIdxVert),2), : ) = [];
% % % 
    
    % set the face color to light gray
    figHandle.FaceColor = [0.75 0.75 0.75];
end