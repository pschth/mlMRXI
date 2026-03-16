function [G, Gx, Gy, Gz] = getVolumetricGradientMatrix(dim)
% computes the volumetric gradient matrices for dimension dim
%%%%
% INPUT:
% dim - 3 x 1 vector; dimensions of the grid to differentiate
%
% OUTPUT:
% G - prod(dim) x prod(dim) matrix; full spatial gradient in all three
% directions (= Gx + Gy + Gz)
% Gx - prod(dim) x prod(dim) matrix; gradient matrix in x-direction
% Gy - prod(dim) x prod(dim) matrix; gradient matrix in y-direction
% Gz - prod(dim) x prod(dim) matrix; gradient matrix in z-direction

%% make Gx
% for gradient in x
Gx = speye(dim(1));
Gx(1:end-1,2:end) = Gx(1:end-1,2:end) - speye(dim(1)-1);
Gx(dim(1),dim(1)) = 0;
% for gradient in y
Gx = kron(eye(dim(2)), Gx);
% for gradient in z
Gx = kron(eye(dim(3)), Gx);

%% make Gy
ny0 = prod(dim(1:2));
ny1 = dim(1)*(dim(2)-1);
Gy = spalloc(ny0, ny0, 2*ny0);
% for gradient in y and x
Gy(1:ny1,:) = spdiags(ones(ny1,1), 0, ny1, ny0) + ...
              spdiags(-ones(ny1,1), dim(1), ny1, ny0);
% for gradient in z
Gy = kron(eye(dim(3)), Gy);

%% make Gz
nz0 = prod(dim);
nz1 = ny0*(dim(3)-1);
Gz = spalloc(prod(dim), prod(dim), 2*prod(dim));
Gz(1:nz1,:) = spdiags(ones(nz1,1), 0, nz1, nz0) + ...
              spdiags(-ones(nz1,1), ny0, nz1, nz0);

%% make G
G = Gx + Gy + Gz;
end

