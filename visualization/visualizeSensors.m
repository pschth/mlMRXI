function sensorHandle = visualizeSensors( SETUP, varargin )
% visualizes the sensor grid
% 
% INPUT:
% SETUP - MRXI setup configuration object
% 
% varargin (optional) - specifies plot properties; has to be stated in
% property pairs (name + value). E.g. varargin{1} = 'Color', varargin{2} =
% 'r'
%
% OUTPUT:
% sensorHandle - 1 x 2 cell; contains the scatter and the quiver plot
% handles for the drawn sensors

% check if properties are given in pairs
if nargin > 1
    nProp = (nargin - 1);
    if mod(nProp,2)~=0
        error('Wrong number of property inputs.');
    end
    nProp = nProp/2;
end

% store initial hold status
holdFlag = ishold;

% plot sensor centers
sensorHandle{1} = scatter3(SETUP.sensorData.centers(1,:),...
                           SETUP.sensorData.centers(2,:),...
                           SETUP.sensorData.centers(3,:),30,'filled');
hold on; % hold for quiver plot
% set varargin properties, if there are any
if nargin > 1
    for i = 1:nProp
        if strcmp(varargin{2*i-1}, 'Color')
            % replace Color property with scatter plot equivalent
            set(sensorHandle{1}, 'MarkerFaceColor', varargin{2*i});
        else
            set(sensorHandle{1}, varargin{2*i-1}, varargin{2*i});
        end
    end
end

% plot sensor orientations
axis tight % shrink displayed figure volume to plotted things
limits = [xlim;ylim;zlim]; % get the current axis limits
quiverScale = norm(diff(limits))/5; % set the quiver arrow scale to a fifth of the axis limit diagonal length
% make quiver plot of sensor orientations
sensorHandle{2} = quiver3(  SETUP.sensorData.centers(1,:),...
                            SETUP.sensorData.centers(2,:),...
                            SETUP.sensorData.centers(3,:),...
                            SETUP.sensorData.orientations(1,:)*quiverScale,...
                            SETUP.sensorData.orientations(2,:)*quiverScale,...
                            SETUP.sensorData.orientations(3,:)*quiverScale);
% set varargin properties, if there are any
if nargin > 1
    for i = 1:nProp
        set(sensorHandle{2}, varargin{2*i-1}, varargin{2*i});
    end
end

% restore initial hold status
if holdFlag
    hold on
else
    hold off;
end
end

