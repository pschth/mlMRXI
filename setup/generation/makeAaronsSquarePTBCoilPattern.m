function coilpattern = makeAaronsSquarePTBCoilPattern( sideLength, nWindings )
% coilpattern = makeAaronsSquarePTBCoilPattern( nWindings )
% 
% generates Aaron's square coil pattern used for the PTB measurements on
% 08/2019 for 1mm boards
% % % 
% INPUT:
% sideLength ... scalar; defines edgelengths of coil
% nWindings ... scalar (optional); defines number of windings. If not
% provided, a wiredistance of approximately 5mil = 0.127mm plus a wirewidth
% of 0.7mm will be set (as in original coils). 
% % % 
% OUTPUT:
% coilpattern ... 3 x nSegments; set of coordinates for the coilpattern.
% nSegments depends on number of windings.

if nargin < 2
    nWindings = sideLength/(2*(0.127e-3+0.7e-3));
elseif isempty(nWindings)
    nWindings = sideLength/(2*(0.127e-3+0.7e-3));
end

% vertical coil distance
zDist = 0.001;

% round nWindings to integers
nWindings = floor(nWindings);

% calculate double wire distance from nWindings to get exact sidelengths
dWireDist = sideLength/(2*(1+nWindings));

% initialize coilpattern
coilpattern = zeros(3,4*nWindings);
coilpattern(:,1) = [-1.5*dWireDist -0.5*dWireDist 0]';

% generate square spiral coil
for w = 1:nWindings
    i = 4*(w-1)+2;
    coilpattern(:,i)    = [-(w+0.5)*dWireDist        (w+0.5)*dWireDist    0]';
    coilpattern(:,i+1)  = [(w+0.5)*dWireDist      (w+0.5)*dWireDist    0]';
    coilpattern(:,i+2)  = [(w+0.5)*dWireDist      -(w+0.5)*dWireDist   0]';
    coilpattern(:,i+3)  = [-(w+1.5)*dWireDist   -(w+0.5)*dWireDist   0]';
end
% offset to origin
offset = coilpattern(:,1) + [dWireDist; dWireDist; 0];
coilpattern = bsxfun(@minus, coilpattern, offset);

% generate bottom part of coil
coilpattern_bottom = [coilpattern(2,:); coilpattern(1,:); coilpattern(3,:)];
coilpattern_bottom = bsxfun(@minus, coilpattern_bottom, coilpattern_bottom(:,end) - coilpattern(:,end));

% unite the two parts with zDist offset
coilpattern = [ bsxfun(@plus, coilpattern(:,end:-1:1), [0 0 zDist/2]') ...
                bsxfun(@minus, coilpattern_bottom, [0 0 zDist/2]')];
end


