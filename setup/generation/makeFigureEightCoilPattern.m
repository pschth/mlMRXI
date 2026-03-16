function [coilpattern, nWindings] = makeFigureEightCoilPattern(...
    innerRadius, outerRadius, solenoidLength, radialWindings, axialWindings, segmentsPerWinding, solenoidCenterDistance)
% [coilpattern, nWindings] = makeSolenoidCoilPattern(...
%     innerRadius, outerRadius, length, radialWindings, axialWindings, nWindings)
% 
% generates a multi-layer double solenoid (figure 8) coil pattern. All
% inputs correspond to ONE (!) of the two solenoids. The other one is
% copied.
% % % 
% INPUT:
% innerRadius ... scalar; defines inner radius of solenoid
% outerRadius ... scalar; defines outer radius of solenoid
% solenoidLength ... scalar; defines length of solenoid
% radialWindings ... scalar; defines number of windings in radial direction
% axialWindings ... scalar; defines number of windings in axial direction
% segmentsPerWinding ... scalar (optional); defines the number of
% discretization segments per winding. Default: 36
% solenoidCenterDistance ... scalar (optional); distance between the two
% outer solenoid borders. Default: 0
% % % 
% OUTPUT:
% coilpattern ... 3 x nSegments; set of coordinates for the coilpattern.
% nSegments depends on number of windings.

if nargin < 6 || isempty(segmentsPerWinding)
    segmentsPerWinding = 36;
end
if nargin < 7 || isempty(solenoidCenterDistance)
    solenoidCenterDistance = 0;
end

% calculate the total number of windings
nWindings = radialWindings*axialWindings;

% calculate total number of segments
nSegments = segmentsPerWinding*nWindings + radialWindings;

% calculate radial, axial and angle increments
if radialWindings == 1
    radialIncrements = innerRadius;
else
    radialIncrements = linspace(innerRadius, outerRadius, radialWindings);
end
axialIncrements = linspace(0, -solenoidLength, segmentsPerWinding*axialWindings+1);
angleIncrements = linspace(0, 2*pi, segmentsPerWinding+1);
angleIncrements = angleIncrements(1:end-1);

% make first solenoid layer
coilpattern = [ radialIncrements(1) .* cos([repmat(angleIncrements, 1, axialWindings) 0]); ...
                radialIncrements(1) .* sin([repmat(angleIncrements, 1, axialWindings) 0]); ...
                axialIncrements];

% make subsequent layers
if radialWindings>1
   unitcoilpattern_odd = [coilpattern(1:2,:)./radialIncrements(1); coilpattern(3,:)];
   unitcoilpattern_even = fliplr(bsxfun(@times, [1 -1 1]', unitcoilpattern_odd));
   
   coilpattern = zeros(3, nSegments);
   segmentsPerLayer = segmentsPerWinding*axialWindings+1;
   for w = 1:radialWindings
       idx = (w-1)*segmentsPerLayer+1:w*segmentsPerLayer;
       if mod(w,2)==1
           % odd
           coilpattern(:,idx) = [unitcoilpattern_odd(1:2,:)*radialIncrements(w); unitcoilpattern_odd(3,:)];
       else
           % even
           coilpattern(:,idx) = [unitcoilpattern_even(1:2,:)*radialIncrements(w); unitcoilpattern_even(3,:)];
       end
   end
   
end

% copy coil pattern to generate second solenoid and center
coilpattern(1,:) = coilpattern(1,:) - (coilpattern(1,end) + solenoidCenterDistance/2);
coilpattern = [coilpattern [-coilpattern(1:2,end:-1:1); coilpattern(3,end:-1:1)]];


