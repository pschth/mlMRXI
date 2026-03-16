function out = scatter3fast(in, varargin)
% uses the scatter3 function on a 3 x X or a X x 3 matrix 'in'
% for a 3 x 3 matrix the x-data has to be in the first row, y-data second
% row, z-data third row

dims = size(in);
if ~any(dims) == 3
    error('Input matrix has no dimension equal to 3.');
end
if dims(1) ~= 3 && dims(2) == 3
    in = in';
end
if nargout > 0
    if nargin > 1
        out = scatter3(in(1,:), in(2,:), in(3,:), varargin{:});
    else
        out = scatter3(in(1,:), in(2,:), in(3,:));
    end
else
    if nargin > 1
        scatter3(in(1,:), in(2,:), in(3,:), varargin{:});
    else
        scatter3(in(1,:), in(2,:), in(3,:));
    end
end
end