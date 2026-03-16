function coilPattern = makeZigzagCoilPattern( coilSize, cornersPerSide, mirror )
% coilPattern = makeZigzagCoilPattern( coilSize, cornersPerSide, (mirror) )
% 
% returns flat, diagonal zig-zag coil pattern.
% 
% % 
% INPUT:
% coilSize - 2 x 1 vector; coil dimensions in x- and y-directions
% cornersPerSide - scalar; defines the corners per coil edge
% mirror - boolean (optional); mirrors the coil along an x-axis. Default:
% false.
% 
% % 
% OUTPUT:
% coilPattern - 3 x (cornersPerSide-1)*4 matrix; standard coil pattern as
% used in the MRXi_IEBE toolbox

%%% exceptions:
if nargin < 3
    mirror = false;
elseif isempty(mirror)
    mirror = false;
end

if length(coilSize) ~= 2
    error('Wrong size for "coilSize". Must be vector with length 2.');
end

if length(cornersPerSide) ~= 1
    error('Wrong size for "cornersPerSide". Must be scalar.');
end

%%% code
% mesh grid for simple corner point extraction
[X,Y] = meshgrid(   linspace(-coilSize(1)/2, coilSize(1)/2, cornersPerSide), ...
                    linspace(-coilSize(2)/2, coilSize(2)/2, cornersPerSide));

% extract corner points for upper and right edges (Northeast) as well as
% for left and lower edges (Southwest) 
stepsNortheast = [X(1,:) X(2:end-1,end)'; Y(1,:) Y(2:end-1,end)'; zeros(1,cornersPerSide*2-2)];
stepsSouthwest = [X(2:end-1,1)' X(end,:); Y(2:end-1,1)' Y(end,:); zeros(1,cornersPerSide*2-2)];

% define matrix column indices for steps
idxNortheast = sort([1:4:cornersPerSide*4-4 2:4:cornersPerSide*4-4]);
idxSouthwest = sort([3:4:cornersPerSide*4-4 4:4:cornersPerSide*4-4]);

% fill coilPattern matrix
coilPattern = zeros(3,cornersPerSide*4-4);
coilPattern(:,idxNortheast) = stepsNortheast;
coilPattern(:,idxSouthwest) = stepsSouthwest;

% mirror coilPattern, if demanded
if mirror
    coilPattern = [coilPattern(1,end:-1:1); coilPattern(2:3,:)];
end
end

