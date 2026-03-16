function [coilpattern, nWindings] = makeSolenoidCoilPattern(...
    innerRadius, outerRadius, solenoidHeight, radialWindings, axialWindings, segmentsPerWinding)
% [coilpattern, nWindings] = makeSolenoidCoilPattern(...
%     innerRadius, outerRadius, length, radialWindings, axialWindings, nWindings)
% 
% generates multi-layer solenoid coil pattern
% % % 
% INPUT:
% innerRadius ... scalar; defines inner radius of solenoid
% outerRadius ... scalar; defines outer radius of solenoid
% solenoidHeight ... scalar; defines length of solenoid
% radialWindings ... scalar; defines number of windings in radial direction
% axialWindings ... scalar; defines number of windings in axial direction
% segmentsPerWinding ... scalar (optional); defines the number of
% discretization segments per winding. Default: 36
% % % 
% OUTPUT:
% coilpattern ... 3 x nSegments; set of coordinates for the coilpattern.
% nSegments depends on number of windings.

if nargin < 6 || isempty(segmentsPerWinding)
    segmentsPerWinding = 36;
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
axialIncrements = linspace(0, -solenoidHeight, segmentsPerWinding*axialWindings+1);
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

