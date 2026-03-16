function coilpattern = makeAntiparallelFFCCoilPattern( R, N2 )
% produces a circular coil pattern with a smaller antiparallel coil within.
% Radius R. It is centered at origin and oriented along the z axis
% % % 
% INPUT:
% R ... scalar, outer radius of the coil
% N2 ... scalar (optional), integer number of windings of inner coil. The
% inner coil radius is defined by the number of windings. Default: 2

if nargin < 2
    N2 = 2;
elseif isempty(N2)
    N2 = 2;
end   

if N2 < 2
    error('N2 cannot be smaller than 2.');
end

N2 = round(N2);

nSegments = 100;
angle = linspace(0,2*pi,nSegments+1);
coilpattern = [ cos(angle);...
                sin(angle);
                zeros(1,length(angle))];
            
coilpattern = R*[coilpattern, repmat([  1/sqrt(N2)*coilpattern(1,:); ...
                                        -1/sqrt(N2)*coilpattern(2,:);...
                                        coilpattern(3,:)],1,N2)];

end

