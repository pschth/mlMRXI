function coilpattern = makeRectangularCoilPattern( nWindings, xDist, yDist, zDist, cw )
% coilpattern = MAKERECTANGULARCOILPATTERN( xDist, yDist, nWindings )
% 
% creates a rectangular coil pattern with orientation [0 0 1], dimensions
% xDist in x- and yDist in y-direction and nWindings windings. Optionally
% this will generate a two layer coil with zDist offset in z-direction. 
% % %
% INPUT:
% nWindings ... scalar; number of windings, can be a real number but will
% be rounded to quarters (e.g. 3.3 -> 3.25)
% xDist ... scalar; distance in x-direction
% yDist ... scalar; distance in y-direction
% zDist ... scalar (optional); distance in z-direction. This will generate
% a two layer coil zDist offset in z-direction.
% cw ... boolean (optional); defines clockwise or counterclockwise spiral
% direction from inside to outside. Default: false
% % % 
% OUTPUT:
% coilpattern ... 3 x nSegments; set of coordinates for the coilpattern.
% nSegments depends on number of windings.

if nargin < 5
    cw = false;
elseif isempty(cw)
    cw = false;
end

% round nWindings to quarters
nWindings = round(nWindings*4)/4;

% half xDist and yDist since the coils are computed from -xDist to +xDist
% (resp. yDist)
xDist = xDist/2;
yDist = yDist/2;

% initialize coilpattern
coilpattern = zeros(3,1+4*nWindings);
coilpattern(:,1) = [xDist/nWindings 0 0]';

% generate rectangular spiral coil
for w = 1:nWindings
    i = 4*(w-1)+2;
    coilpattern(:,i)    = [w*xDist/nWindings        w*yDist/nWindings    0]';
    coilpattern(:,i+1)  = [-w*xDist/nWindings       w*yDist/nWindings    0]';
    coilpattern(:,i+2)  = [-w*xDist/nWindings       -w*yDist/nWindings   0]';
    coilpattern(:,i+3)  = [(w+1)*xDist/nWindings    -w*yDist/nWindings   0]';
end
% continue spiral for quarter turns
if nWindings - floor(nWindings) > eps
    w = w+1;
    i = 4*(w-1)+2;
    coilpattern(:,i)    = [w*xDist/nWindings        w*yDist/nWindings    0]';
    
    if nWindings - floor(nWindings) - 1/4 > eps
        coilpattern(:,i+1)  = [-w*xDist/nWindings       w*yDist/nWindings    0]';
        
        if nWindings - floor(nWindings) - 1/2 > eps
            coilpattern(:,i+2)  = [-w*xDist/nWindings       -w*yDist/nWindings   0]';
        end
    end
end

% mirror coil to clockwise rotation if demanded
if cw
    coilpattern(1,:) = -coilpattern(1,:);
end


% if zDist is provided
if nargin > 3
    if ~isempty(zDist)
        if zDist ~= 0
            % initialize bottom coilpattern
            coilpattern_bottom = zeros(3,1+4*nWindings);
            coilpattern_bottom(:,1) = [xDist/nWindings 0 0]';
            coilpattern_bottom(:,2) = [-xDist/nWindings 0 0]';

            % generate rectangular bottom spiral coil
            for w = 1:nWindings
                i = 4*(w-1)+3;
                coilpattern_bottom(:,i)    = [-w*xDist/nWindings        w*yDist/nWindings    0]';
                coilpattern_bottom(:,i+1)  = [(w+1)*xDist/nWindings       w*yDist/nWindings    0]';
                coilpattern_bottom(:,i+2)  = [(w+1)*xDist/nWindings       -w*yDist/nWindings   0]';
                coilpattern_bottom(:,i+3)  = [-(w+1)*xDist/nWindings    -w*yDist/nWindings   0]';
            end
            
            % continue bottom spiral for quarter turns until it terminates
            % beneath the upper coil beginning
            if nWindings - floor(nWindings) < eps
                
                coilpattern_bottom(:,end) = [];
                
            else
                w = w+1;
                i = 4*(w-1)+3;
                
                if nWindings - floor(nWindings) > eps
                    coilpattern_bottom(:,i)  = [-w*xDist/nWindings       w*yDist/nWindings    0]';
                    coilpattern_bottom(:,i+1)  = [w*xDist/nWindings       w*yDist/nWindings    0]';
                    if nWindings - floor(nWindings) - 1/4 > eps
                        coilpattern_bottom(:,i+1)  = [];
                        if nWindings - floor(nWindings) - 1/2 > eps
                            coilpattern_bottom(:,i+1)  = [(w+1)*xDist/nWindings       w*yDist/nWindings    0]';
                            coilpattern_bottom(:,i+2)  = [(w+1)*xDist/nWindings       -w*yDist/nWindings    0]';
                            coilpattern_bottom(:,i+3)  = [-w*xDist/nWindings       -w*yDist/nWindings    0]';
                        end
                    end
                end
            end

            % mirror coil to clockwise rotation if demanded
            if cw
                coilpattern_bottom(1,:) = -coilpattern_bottom(1,:);
            end
            
            % unite the two parts with zDist offset
            coilpattern = [ bsxfun(@plus, coilpattern(:,end:-1:1), [0 0 zDist/2]') ...
                            bsxfun(@minus, coilpattern_bottom, [0 0 zDist/2]')];
        end
    end
end


end

