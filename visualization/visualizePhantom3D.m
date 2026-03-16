function handles = visualizePhantom3D( mrxisetup, phantoms, useCommonCRange, numberOfColumns )
% visualizes an MNP distribution defined in the vector 'phantom' based on
% the ROIData defined in 'mrxisetup'. 
% ATTENTION: works only if mrxisetup.ROIData.voxelDims is defined.
% 
% INPUT:
% mrxisetup - standard mrxi SETUP object
% 
% phantoms - nVoxel x 1 vector or cell of nVoxel x 1 vectors; contains the
% MNP concentrations per voxel. Multiple MNP concentrations can be passed
% in a cell.
%
% useCommonCRange - boolean (optional); if TRUE, a common color range is
% used among all phantoms. Default: FALSE
% 
% numberOfColumns - scalar (optional); specifies number of subfigure
% columns in output figure. Default: min(numberOfPassedPhantoms, 3)
% 
% OUTPUT:
% handles - 3 x numberOfPassedPhantoms plot handles; 
% 1st row: patch handle(s) for visualizing the outlines of the voxel grid
% 2nd row: patch handle(s) for visualizing the MNP distribution
% 3rd row: axis handle(s) for the subfigure axis/axes
% 
% 
% HINT: the axis link property is stored in the figure properties under
% gcf.UserData.AxisLink

if nargin < 3 || isempty(useCommonCRange)
    useCommonCRange = false;
end

% check if phantoms is cell, otherwise convert to cell for consistent
% subsequent computations
if ~iscell(phantoms)
    temp = phantoms;
    phantoms = cell(1);
    phantoms{1} = temp;
end
% get number of phantoms to visualize
nPhants = length(phantoms);

% check if numberOfColumns is given, otherwise pass default value
if nargin < 4 || isempty(numberOfColumns)
    numberOfColumns = min(nPhants,3);
end

% get hold status and clear figure if hold is 'off'
holdstatus = ishold;
if ~holdstatus
    gcf; clf;
end

% get the minimum and maximum voxel corner points of a voxel centered at
% origin
l = -mrxisetup.ROIData.voxelDims./2;
u = mrxisetup.ROIData.voxelDims./2;

% define patch coordinates centered at origin: per coordinate 6 rows (6
% faces of the cuboidal voxel) and 4 rows (4 corners per face)
X = [l(1) l(1) l(1) l(1); l(1) l(1) u(1) u(1); l(1) l(1) u(1) u(1);...
     u(1) u(1) u(1) u(1); l(1) l(1) u(1) u(1); l(1) l(1) u(1) u(1)];
Y = [l(2) l(2) u(2) u(2); l(2) l(2) l(2) l(2); l(2) u(2) u(2) l(2);...
     l(2) l(2) u(2) u(2); u(2) u(2) u(2) u(2); l(2) u(2) u(2) l(2)];
Z = [l(3) u(3) u(3) l(3); l(3) u(3) u(3) l(3); l(3) l(3) l(3) l(3);...
     l(3) u(3) u(3) l(3); l(3) u(3) u(3) l(3); u(3) u(3) u(3) u(3)];

% rotate patch coordinates to orientation of the ROI
pCoords = rotatepoints([X(:) Y(:) Z(:)], mrxisetup.ROIData.rotvector);

% get number of voxels
nVoxels = mrxisetup.getNumberOfVoxels;
% initialize patch matrices for all voxels for x-, y-, and z-coordinates
XAll = zeros(4, 6*nVoxels);
YAll = XAll;
ZAll = XAll;
% for each voxel...
for v = 1:nVoxels
    % calculate fill index for the patch matrices
    fillIdx = (v-1)*6+1:v*6;
    % translate rotated patch coordinates (still centered at origin) to
    % actual voxel positions
    pCoordsCurr = bsxfun(@plus, pCoords, mrxisetup.ROIData.voxels(:,v));
    % store rotated and translated patch coordinates for every voxel
    XAll(:, fillIdx) = reshape(pCoordsCurr(1,:), 6, 4)';
    YAll(:, fillIdx) = reshape(pCoordsCurr(2,:), 6, 4)';
    ZAll(:, fillIdx) = reshape(pCoordsCurr(3,:), 6, 4)';
end

% calculate number of rows of subfigures
numberOfRows = ceil(nPhants/numberOfColumns);
% initialize handle storage cells
phantHandle = cell(1, nPhants);
roiHandle = phantHandle;
axHandle = phantHandle;
% if common color range is required, the minimum and maximum color range of
% all phantoms has to be extracted. For that, initialize the color limits
% with cRange = [lowerLimit=Inf, upperLimit=-Inf] and successively update
% with every visualized phantom
if useCommonCRange
    cRange = [Inf -Inf];
end
% for every phantom...
for p = 1:length(phantoms)
    % make a new subfigure in figure
    % INFO FOR pyMRXI: I guess something like 'subaxis' should be possible
    % by using matplotlib...
    axHandle{p} = subaxis(numberOfRows, numberOfColumns, p, ...
        'Spacing',0.0, 'Padding',0.03,'Margin',0.0);
    
    if p < 2
        % draw outline of ROI using 'drawvoxels' for first phantom
        roiHandle{p} = drawvoxels(mrxisetup);
        % change some visualization settings to make the outline less
        % visible
        roiHandle{p}.FaceAlpha = 0;
        roiHandle{p}.LineStyle = ':';
        roiHandle{p}.EdgeColor = [1 1 1]*0.75;
        roiHandle{p}.EdgeAlpha = 0.5;
    else
        % copy outline of ROI from first phantom to subsequent phantoms to
        % avoid computing it again
        roiHandle{p} = copyobj(roiHandle{1}, axHandle{p});
    end
    
    % only visualize voxels with an MNP mass not equal to 0. This boolean
    % vector is reproduced 6 times for the 6 faces of a single voxel.
    useVoxels = repmat(phantoms{p}(:)' ~= 0, 6, 1);
    % bring into column vector format for indexing
    useVoxels = useVoxels(:);
    % do the same thing for the color data (=MNP mass per voxel)
    CData = repmat(phantoms{p}(:)', 6, 1);
    CData = CData(:);
    
    % sort according to CData (opacity at adjacent voxels is plot order
    % dependent and I want voxels with more MNP content to be more visible
    % than voxels with fewer MNPs)
    patchData = sortrows([XAll(:,useVoxels)' YAll(:,useVoxels)' ZAll(:,useVoxels)' CData(useVoxels)], 13);
    
    % visualize voxels as patches
    phantHandle{p} = patch( patchData(:,1:4)', ... voxel patch x-coordinates
                            patchData(:,5:8)', ... voxel patch y-coordinates
                            patchData(:,9:12)', ... voxel patch z-coordinates
                            patchData(:,13)); % voxel patch color data
    
    % get minimum color data of all voxels containing MNPs
    minCData = min(CData(useVoxels));
    % use zero as minimum if CData(useVoxels) contains only positive values
    minCData = min(minCData,0);
    % get maximum color data of all voxels containing MNPs
    maxCData = max(CData(useVoxels));
    % if minimum and maximum color data are not equal...
    if minCData ~= maxCData
        % set color scale and opacity scale to range from minimum to
        % maximum color data values
        axHandle{p}.CLim = [minCData maxCData];
        axHandle{p}.ALim = [minCData maxCData];
        % add opacity data to the phantom visualization (voxels with less
        % MNP are more see-through than voxels with more MNPs)
        phantHandle{p}.FaceVertexAlphaData = patchData(:,13);
        % make opacity map from barely visible to complete opacity in 64
        % steps
        phantHandle{p}.Parent.Alphamap = linspace(0.05, 1, 64);
        % scale actual color data to provided Alphamap
        phantHandle{p}.FaceAlpha = 'flat';
        phantHandle{p}.EdgeAlpha = 'flat';
        
    % if minimum and maximum color data are equal...
    else
        % set color scale and opacity scale to a little bit below maximum
        % and maximum color data (this way, only the maxCData color and
        % opacity will be displayed for the MNP filled voxels)
        axHandle{p}.CLim = [maxCData-eps maxCData];
        % set MNP filled voxels to full opacity
        phantHandle{p}.FaceAlpha = 1;
        phantHandle{p}.EdgeAlpha = 1;
    end
 
    % set all axis scales to one to avoid image distortion
    daspect([1 1 1]);
    % set to a typically good view of the phantoms
    view([30 40]);
    % label axes
    xlabel x; ylabel y; zlabel z;
    % set axis limits to tightly fit to subfigure content
    axis tight;
    
    % if a common color range is required...
    if useCommonCRange
        % get the minimum and maximum values of the stored color range and
        % the currently plotted color range
        cRange = [min(cRange(1), axHandle{p}.CLim(1))  max(cRange(2), axHandle{p}.CLim(2))];
        % use zero as minimum if axHandle{p}.CLim(1) contains only positive
        % values
        cRange(1) = min(cRange(1), 0);
    end
    % store axis handle
    axVec(p) = axHandle{p}; %#ok<AGROW>
end

% if common color range is required, set the color and opacity limits of
% all axes to the same color range
if useCommonCRange
    for p = 1:length(phantoms)
        axHandle{p}.CLim = cRange;
        axHandle{p}.ALim = cRange;
    end
end

% link camera positions of all subfigures to be able to rotate all
% subfigures simultaneously to the same degree
f = gcf;
f.UserData.AxisLink = linkprop(axVec,{'CameraPosition','CameraUpVector'});

% define color map (ranges from white over bright red to dark red)
cmap = [[linspace(1, 1, 256); linspace(1, 0, 256); linspace(1, 0, 256)] ...
       [linspace(1, 0.75, 128); linspace(0, 0, 128); linspace(0, 0, 128)]]';
colormap(cmap);

% restore initial hold status
if holdstatus
    hold on
else
    hold off
end

% return handles
handles = [roiHandle; phantHandle; axHandle];
end

