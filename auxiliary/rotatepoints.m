function [newpoints, rotMatrix] = rotatepoints(points, direction, startVec)
% Rotates a pointcloud "points" with initial orientation along the z-axis
% centered at the origin to a new directional vector "direction". If
% startVec is given, the rotation is not initialized along the z-axis, but
% along startVec.
%
% INPUT:
% points - 3 x nPoints matrix; pointcloud to be rotated around the origin
% and originally oriented along the z-direction
%
% direction - 3 x 1 vector; points along the new orientation of the
% pointcloud
% 
% startVec - 3 x 1 vector (optional); initial orientation of the points,
% Default: [0 0 1]'
%     
% OUTPUT: newpoints - 3 x nPoints matrix; from startVec to direction
% rotated pointcloud
% 
% rotMatrix - 3 x 3 matrix; rotation matrix so that newpoints = rotMatrix *
% points

    %%% exceptions
    % check input sizes
    if ~any(size(points) == 3)
        error('Wrong size of "points" matrix.');
    end
    if ~(numel(direction) == 3)
        error('Wrong size of "direction" vector.');
    end
    % check if startVec is given
    if nargin < 3 || isempty(startVec)
        startVec = [0 0 1]';
    else
        % normalize startVec vector, if given
        startVec = normcols(startVec(:));
    end

    % transpose pointcloud matrix if necessary
    if size(points,1) ~= 3 && size(points,2) == 3
        points = points';
    end
    
    % normalize direction vector
    direction = normcols(direction(:));

    %%% code
    % exception: direction equals z-vector or negative z-vector
    if all(direction==startVec)
        rotMatrix = eye(3);
    elseif all(direction==-startVec)
        rotMatrix = -eye(3);  
    else
        % get rotation axis (normalized cross product of startVec with direction)
        rotAxis = normcols(cross(startVec,direction));
        % get rotation angle
        rotAngle = acos(dot(startVec,direction));
        % calculate rotation matrix
        rotMatrix = [  rotAxis(1)*rotAxis(1)*(1-cos(rotAngle))+cos(rotAngle)                   rotAxis(1)*rotAxis(2)*(1-cos(rotAngle))-rotAxis(3)*sin(rotAngle)        rotAxis(1)*rotAxis(3)*(1-cos(rotAngle))+rotAxis(2)*sin(rotAngle);
                       rotAxis(2)*rotAxis(1)*(1-cos(rotAngle))+rotAxis(3)*sin(rotAngle)        rotAxis(2)*rotAxis(2)*(1-cos(rotAngle))+cos(rotAngle)                   rotAxis(2)*rotAxis(3)*(1-cos(rotAngle))-rotAxis(1)*sin(rotAngle);
                       rotAxis(3)*rotAxis(1)*(1-cos(rotAngle))-rotAxis(2)*sin(rotAngle)        rotAxis(3)*rotAxis(2)*(1-cos(rotAngle))+rotAxis(1)*sin(rotAngle)        rotAxis(3)*rotAxis(3)*(1-cos(rotAngle))+cos(rotAngle) ];
    end
    % rotate pointcloud
    newpoints = rotMatrix*points; 
end