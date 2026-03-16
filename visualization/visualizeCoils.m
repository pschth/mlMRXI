function coilHandle = visualizeCoils( mrxisetup, varargin )
% visualizes the coil grid of a SETUP object
% 
% INPUT:
% mrxisetup - MRXI setup configuration object
% 
% varargin (optional) - specifies plot properties; has to be stated in
% property pairs (name + value). E.g. varargin{1} = 'Color', varargin{2} =
% 'r'
%
% OUTPUT:
% coilHandle - nCoils x 1 cell; contains the plot handles for all drawn
% coils

% check if object has coil centers, otherwise return empty handle
if isempty(mrxisetup.coilData.centers)
    coilHandle = [];
    return
end

% check if properties are given in pairs
if nargin > 1
    nProp = (nargin - 1); % get number of elements in 'varargin'
    if mod(nProp,2)~=0 % check if number of elements is even
        error('Wrong number of property inputs.');
    end
    nProp = nProp/2; % get number of property pairs
end

% store hold status
holdFlag = ishold;

% if coils are given as dipoles, visualize them as scatter and quiver plots
if mrxisetup.coilData.isDipole
    % scatter plot of coil centers
    coilHandle{1} = scatter3(mrxisetup.coilData.centers(1,:),mrxisetup.coilData.centers(2,:),mrxisetup.coilData.centers(3,:),30,'filled','r');
    hold on;
    % set properties, if any are given
    if nargin > 1
        for i = 1:nProp
            if strcmp(varargin{2*i-1}, 'Color')
                % replace Color property with scatter plot equivalent
                set(coilHandle{1}, 'MarkerFaceColor', varargin{2*i});
            else
                set(coilHandle{1}, varargin{2*i-1}, varargin{2*i});
            end
        end
    end
    
    % quiver plot for coil orientations
    if size(mrxisetup.coilData.centers,2)==1
        % if only one coil is being visualized calculate maximum vector
        % length (half of ROI diagonal)
        l = norm(mrxisetup.ROIData.ROI(:,2) - mrxisetup.ROIData.ROI(:,1))/2; % half the length of the ROI diagonal
        % make quiver plot with calculated arrow lengths
        coilHandle{2} = quiver3(    mrxisetup.coilData.centers(1,:),mrxisetup.coilData.centers(2,:),mrxisetup.coilData.centers(3,:),...
                                    mrxisetup.coilData.orientations(1,:),mrxisetup.coilData.orientations(2,:),mrxisetup.coilData.orientations(3,:),...
                                    l,'AutoScale','off','Color','r');
    else
        % else make quiver plot with autoscale
        coilHandle{2} = quiver3(    mrxisetup.coilData.centers(1,:),mrxisetup.coilData.centers(2,:),mrxisetup.coilData.centers(3,:),...
                                    mrxisetup.coilData.orientations(1,:),mrxisetup.coilData.orientations(2,:),mrxisetup.coilData.orientations(3,:),...
                                    'Color','r');
    end
    % set property pairs to quiver plot
    if nargin > 1
        for i = 1:nProp
            set(coilHandle{2}, varargin{2*i-1}, varargin{2*i});
        end
    end
    
% if coils are no dipoles, plot coil segments    
else
    % get number of coils to plot
    nCoils = mrxisetup.getNumberOfCoils;
    
    % generate different shades of red for different coil patterns from hsv
    % interpolation if there are more than one coil pattern defined
    nCoilpatterns = length(mrxisetup.coilData.coilpattern);
    if nCoilpatterns > 1
        hInterpVals = [0 0]; % hue: keep at zero (=red)
        sInterpVals = [1 1 0.3]; % saturation: keep at one for first half of coil patterns, then go linearly down to 0.3 for the second half
        vInterpVals = [0.2 1 1]; % value: go linearly from 0.2 to one for the first half of coil patterns, then keep at one for the second half
        % compute the gradients like described above
        hGradient = interp1(linspace(0,1,length(hInterpVals)), hInterpVals, linspace(0,1,nCoilpatterns));
        sGradient = interp1(linspace(0,1,length(sInterpVals)), sInterpVals, linspace(0,1,nCoilpatterns));
        vGradient = interp1(linspace(0,1,length(vInterpVals)), vInterpVals, linspace(0,1,nCoilpatterns));
        % merge into matrix of HSV color values
        hsvVals = [hGradient; sGradient; vGradient]';
        % convert from HSV to RGB
        rgbVals = hsv2rgb(hsvVals);
    else
        % use red coils if only one coil pattern is used
        rgbVals = [1 0 0];
    end
    
% % % % % % % to catch exceptions from old SETUP objects. INFO FOR pyMRXI: not
% % % % % % % necessary!
    if ~isfield(mrxisetup.coilData, 'coilpatternassignment')
        mrxisetup.coilData.coilpatternassignment = ones(mrxisetup.getNumberOfCoils,1);
    end
% % % % % % %     

    % initialize plot3 handle storage cell
    coilHandle = cell(nCoils,1);
    % loop through all coils
    
%             %%%%%%%% TEMPORARY: JUST FOR JOURNAL PAPER FIGURE
%             rgbVals = colormap('parula');
%             idx = round(linspace(1,size(rgbVals,1), nCoilpatterns));
%             rgbVals = rgbVals(idx,:);
%             bwVals8 = repmat([0.1 0.1 0.1; 0.7 0.7 0.7],4,1);
%             rgbVals = repmat([bwVals8; circshift(bwVals8,1,1)],2,1);
%             alphaArray = linspace(1,1,nCoils);
%             %%%%%%%%
            
    for c = 1:nCoils
        % get name of current filamentary segment field in mrxisetup.coilData.segments
        segName = ['c' num2str(c)];
        lw = 0.5; % plot line width
        

        
        % 3D plot of coils with defined color of current coilpattern and
        % the defined line width
        coilHandle{c} = plot3(  mrxisetup.coilData.segments.(segName)(1,:),...
                                mrxisetup.coilData.segments.(segName)(2,:),...
                                mrxisetup.coilData.segments.(segName)(3,:),...
                                'Color', rgbVals(mrxisetup.coilData.coilpatternassignment(c),:),...
                                'LineWidth', lw);
                            
%             %%%%%%%% TEMPORARY: JUST FOR JOURNAL PAPER FIGURE
%             [k1,~] = convhull(mrxisetup.coilData.segments.(segName)(1,:),...
%                                 mrxisetup.coilData.segments.(segName)(2,:),...
%                                 mrxisetup.coilData.segments.(segName)(3,:));
%             coilHandle{c} = ...
%                 trisurf(k1, mrxisetup.coilData.segments.(segName)(1,:),...
%                             mrxisetup.coilData.segments.(segName)(2,:),...
%                             mrxisetup.coilData.segments.(segName)(3,:),...
%                             'FaceColor',rgbVals(mrxisetup.coilData.coilpatternassignment(c),:),...
%                             'EdgeColor',rgbVals(mrxisetup.coilData.coilpatternassignment(c),:),...
%                             'FaceAlpha',alphaArray(c), 'EdgeAlpha', 0);
%             
%             %%%%%%%%
        hold on; % hold figure for plotting multiple coils
        % set varargin properties, if there are any
        if nargin > 1
            for i = 1:nProp
                set(coilHandle{c}, varargin{2*i-1}, varargin{2*i});
            end
        end
    end
end
% return to initial hold status
if holdFlag
    hold on;
else
    hold off;
end
end

