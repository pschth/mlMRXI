function [coilpattern, nWindings] = makeSpiralCoilPattern(radius, nWindings, twoLayerDistance)
% coilpattern = makeSpiralCoilPattern(radius, nWindings)
% 
% generates spiral coil pattern similar to Maik's spiral coil pattern
% % % 
% INPUT:
% radius ... scalar; defines radius of coil start point
% nWindings ... scalar (optional); defines number of windings. If not
% provided, a wiredistance of approximately 5mil = 0.127mm plus a wirewidth
% of 0.7mm will be set (as in original coils). 
% twoLayerDistance ... scalar (optional); distance between bottom and top
% layer of coil. If twoLayerDistance==0 a single layer coil is generated.
% Default twoLayerDistance=0.
% % % 
% OUTPUT:
% coilpattern ... 3 x nSegments; set of coordinates for the coilpattern.
% nSegments depends on number of windings.

if nargin < 2
    nWindings = round(radius/(0.127e-3+0.7e-3));
elseif isempty(nWindings)
    nWindings = round(radius/(0.127e-3+0.7e-3));
end

if nargin < 3
    twoLayerDistance = 0;
elseif isempty(nWindings)
    twoLayerDistance = 0;
end

segmentsPerWinding = 50;
nSegments = segmentsPerWinding*nWindings;
angleInc = linspace(0, nWindings*2*pi, nSegments);
rInc = linspace(0, radius, nSegments);

% make flat spiral pattern
coilpattern = [rInc .* cos(angleInc); rInc .* -sin(angleInc)];
% continue pattern on other side of circuit board
coilpatternBottom = [-coilpattern(1,:); coilpattern(2,:)];
coilpattern = [ [coilpattern(:,end:-1:1); zeros(1,nSegments)] ...
                [0;0;-twoLayerDistance] ...
                [coilpatternBottom; -ones(1,nSegments)*twoLayerDistance] ] + [0;0;twoLayerDistance/2];
% % coilpattern polygon has to be closed for correct computation of the
% % magnetic field
% coilpattern(:,end+1) = coilpattern(:,1);