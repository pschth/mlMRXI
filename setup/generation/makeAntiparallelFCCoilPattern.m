function [coilPattern, rn, fval] = makeAntiparallelFCCoilPattern( R, N1, N2, x, R2init)
% generates a coil pattern with two connected coils that compensates the
% magnetic field at distance x
% % % 
% INPUT:
% R ... scalar, radius of coil 1
% N1, N2 ... scalars, integer number of coil turns of coils 1, 2, and 3.
% Can also be negative to invert coil winding direction.
% x... scalar, distances to field-free points along the central rotatory
% coil axes
% R2init ... scalar, initial guess for coil radius of compensation coil
% % % 
% OUTPUT:
% coilPattern ... 3 x nSegments matrix, calculated coil pattern
% rn ... 2 x 2 matrix, 1st column = coil radius scaling factors, 2nd column
% = respective number of coil turns
% fval ... final function value of minimization to determine coil radii.
% Must be close to zero or minimization most likely failed.
% % %
% VALUES TO TRY:


% define function to null magnetic field at distances x1 and x2
rfun = @(R2)     (N1*R^2/sqrt(R^2+x^2)^3 + N2*R2^2/sqrt(R2^2+x^2)^3)^2;

% find coil radii for given function
opts.TolFun = eps;
[R2, fval] = fminsearch(rfun, R2init, opts);

% warn if function does not converge as intended
if fval > eps
    warning('Function value did not converge below tolerance.');
end

% sort radii from smallest to largest radius
rn = [R N1; R2 N2];
rn = sortrows(rn,1);

% generate a single coil turn with 100 coil segments
nSegments = 100;
angle = linspace(0,2*pi,nSegments+1);
coilTurn = [ cos(angle);...
                sin(angle);
                zeros(1,length(angle))];

% combine single coil turns to produce compensation coil
coilPattern =   [ rn(2,1)*repmat([coilTurn(1,:); sign(rn(2,2))*coilTurn(2,:); coilTurn(3,:)],1,abs(rn(2,2))) ...
                  rn(1,1)*repmat([coilTurn(1,:); sign(rn(1,2))*coilTurn(2,:); coilTurn(3,:)],1,abs(rn(1,2))) ...
                  rn(2,1)*coilTurn(:,1)];
end

