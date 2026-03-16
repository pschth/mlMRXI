function out = quiver3fast(pos, dir, varargin)
% uses the quiver3 function on two 3 x X or X x 3 matrices 'pos' and 'dir'.
% for a 3 x 3 matrix the x-data has to be in the first row, y-data second
% row, z-data third row

dims = size(pos);
if ~any(dims) == 3
    error('Input matrix has no dimension equal to 3.');
end
if dims(1) ~= 3 && dims(2) == 3
    pos = pos';
end
dims = size(dir);
if ~any(dims) == 3
    error('Input matrix has no dimension equal to 3.');
end
if dims(1) ~= 3 && dims(2) == 3
    dir = dir';
end

if nargout > 0
    if nargin > 1
        out = quiver3(pos(1,:), pos(2,:), pos(3,:), ...
                      dir(1,:), dir(2,:), dir(3,:), varargin{:});
    else
        out = quiver3(pos(1,:), pos(2,:), pos(3,:), ...
                      dir(1,:), dir(2,:), dir(3,:));
    end
else
    if nargin > 1
        quiver3(pos(1,:), pos(2,:), pos(3,:), ...
                dir(1,:), dir(2,:), dir(3,:), varargin{:});
    else
        quiver3(pos(1,:), pos(2,:), pos(3,:), ...
                dir(1,:), dir(2,:), dir(3,:));
    end
end
end