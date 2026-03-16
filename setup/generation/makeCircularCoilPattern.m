function coilpattern = makeCircularCoilPattern( R, nSegments )
% produces a circular coil pattern with nSegments coil segments and radius
% R. It is centered at origin and oriented along the z axis
% 
% INPUT:
% R - scalar; radius of coil in m
% 
% nSegments - scalar; number of segments the coil is discretized
% into. Must be > 2.
% 
% OUTPUT:
% coilpattern - 3 x (nSegments+1) matrix; contains the end points of the
% filamentary segments of the circular coil discretization

if nSegments < 3
    % check if enough segments are given to form a circular polygon
    error('Cannot construct 3D structure with this number of nSegments. Has to be > 2.');
else
    % calculate angular increments based on desired number of coil segments
    angle = linspace(0,2*pi,nSegments+1);
    
    % calculate coil segment end point coordinates for the coilpattern
    coilpattern = [ R*cos(angle);...
                    R*sin(angle);
                    zeros(1,length(angle))];
end
end

