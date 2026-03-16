function coilpattern = makeAaronsRectangularCoilPattern( withConnectors )
% coilpattern = makeAaronsRectangularCoilPattern( withConnectors )
% 
% generates the coil pattern used by Aaron
% % % 
% INPUT:
% withConnectors ... boolean (optional); decides if the connector lines are
% in the coil pattern or not. Default: false
% % % 
% OUTPUT:
% coilpattern ... 3 x nSegments; set of coordinates for the coilpattern.
% nSegments depends on number of windings.

if nargin == 0
    withConnectors = false;
elseif isempty(withConnectors)
    withConnectors = false;
end

% generate rough coil pattern 
nWindings = 19.5;
xDist = 0.052;
yDist = 0.052;
zDist = 0.0011;
cw = false;

coilpattern = makeRectangularCoilPattern( nWindings, xDist, yDist, zDist, cw );

% refine coil pattern to exact design
coilpattern(:,1:4) = [];
coilpattern(1,1) = coilpattern(1,1) - xDist/floor(nWindings);
if withConnectors
    coilpattern = [coilpattern(:,1)+[0 2*yDist/floor(nWindings) 0]' coilpattern];
    coilpattern(2,end) = coilpattern(2,end) + yDist/floor(nWindings);
else
    coilpattern(2,end) = coilpattern(2,end) - yDist/floor(nWindings);
end


end

