function rotMatrix = getRotationMatrix(roll, pitch, yaw)
% returns a 3 x 3 rotation matrix defined by yaw, pitch, and roll given in
% [rad]
rotMatrix = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1] * ...
            [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)] * ...
            [1 0 0; 0 cos(roll) -sin(roll); 0 sin(roll) cos(roll)];
end

