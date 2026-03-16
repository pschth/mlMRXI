function [coilPattern, rn, fval] = makeNearAndFarFieldCompensatingCoil( R, N1, N2, N3, x1, x2, r2init, r3init)
% generates a coil pattern with three connected cylindrical coils that
% nulls the magnetic field at distances x1 and x2 by choosing two of the
% three coil radii with a minimization approach. The first coil is fixed at
% radius R with N1 coil turns. 
% % % 
% INPUT:
% R ... scalar, radius of coil 1
% N1, N2, N3 ... scalars, integer number of coil turns of coils 1, 2, and
% 3. Can also be negative to invert coil winding direction.
% x1, x2 ... scalars, distances to field-free points along the central rotatory coil
% axes 
% r2init, r3init ... scalars, initial guess for coil radius scaling
% factors, so that the radii of coils 2 and 3 are equal to r2init*R and
% r3init*R
% % % 
% OUTPUT:
% coilPattern ... 3 x nSegments matrix, calculated coil pattern
% rn ... 3 x 2 matrix, 1st column = coil radius scaling factors, 2nd column
% = respective number of coil turns
% fval ... final function value of minimization to determine coil radii.
% Must be close to zero or minimization most likely failed.
% % %
% VALUES TO TRY:
% R = 0.01;
% N1 = -2;
% N2 = 2;
% N3 = 1;
% x1 = 0.01;
% x2 = 0.13;
% r2init = 0.5;
% r3init = 2;


% define function to null magnetic field at distances x1 and x2
rfun = @(r)     (N1/sqrt(R^2+x1^2)^3 + ...
                N2*r(1)^2/sqrt(r(1)^2*R^2+x1^2)^3 + ...
                N3*r(2)^2/sqrt(r(2)^2*R^2+x1^2)^3)^2 + ...
                (N1/sqrt(R^2+x2^2)^3 + ...
                N2*r(1)^2/sqrt(r(1)^2*R^2+x2^2)^3 + ...
                N3*r(2)^2/sqrt(r(2)^2*R^2+x2^2)^3)^2;

% find coil radii for given function
opts.TolFun = eps;
opts.MaxFunEvals = 1e6;
[r, fval] = fminsearch(rfun, [r2init r3init], opts);

% warn if function does not converge as intended
if fval > eps
    warning('Function value did not converge below tolerance.');
end

% sort radii from smallest to largest radius
rn = [r'*R [N2;N3]; R N1];
rn = sortrows(rn,1);

% generate a single coil turn with 100 coil segments
nSegments = 100;
angle = linspace(0,2*pi,nSegments+1);
coilTurn = [ cos(angle);...
                sin(angle);
                zeros(1,length(angle))];

% combine single coil turns to produce compensation coil
coilPattern =   [ rn(3,1)*repmat([coilTurn(1,:); sign(rn(3,2))*coilTurn(2,:); coilTurn(3,:)],1,abs(rn(3,2))) ...
                  rn(2,1)*repmat([coilTurn(1,:); sign(rn(2,2))*coilTurn(2,:); coilTurn(3,:)],1,abs(rn(2,2))) ...
                  rn(1,1)*repmat([coilTurn(1,:); sign(rn(1,2))*coilTurn(2,:); coilTurn(3,:)],1,abs(rn(1,2))) ...
                  rn(3,1)*coilTurn(:,1)];
end

