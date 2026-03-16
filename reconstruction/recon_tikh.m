function result = recon_tikh( L, b, alpha, Gamma, inv_LtL_alphaxgamma )
% Tikhonov reconstruction
% 
% INPUT:
% L - nSensors*nActivations x nVoxel matrix; leadfield matrix
% 
% b - nSensors*nActivations vector; measurement vector
% 
% alpha - scalar; regularization parameter
% 
% Gamma (optional) - nVoxel vector or nVoxel x nVoxel matrix; regularization matrix
% 
% inv_LtL_alphaxgamma (optional) - nVoxel x nVoxel matrix; inverted term of the
% tikhonov reconstruction: inv(L'*L + alpha*eye(nVoxel)). Can be used for
% faster calculation of multiple reconstructions; default: []

% check L and b
if size(L,1) ~= size(b,1)
    error('L or b has wrong size. (L - nSensors*nActivations x nVoxel matrix, b - nSensors*nActivations vector)');
end

nVox = size(L,2);

if nargin == 3
    % if L, b and alpha are provided: Gamma = eye(nVox)
    
    Gamma = speye(nVox);
    result = (L'*L + alpha^2*Gamma)\L'*b;
    
elseif nargin == 4 && ~isempty(Gamma)
    % if Gamma is provided additionally
    
    % check right size of Gamma
    if ~(all(size(Gamma)==nVox) || all(size(Gamma)==[nVox 1]) || all(size(Gamma)==[1 nVox]))
        error('Gamma has wrong size. (Gamma - nVoxel vector or nVoxel x nVoxel matrix)');
    elseif all(size(Gamma)==[nVox 1]) || all(size(Gamma)==[1 nVox])
        Gamma = sparse(1:nVox, 1:nVox, Gamma(:));
    end
    
    result = (L'*L + alpha^2*Gamma)\(L'*b);
    
elseif nargin == 5
    % if inv_LtL_alphaxgamma is provided additionally
    
    result = inv_LtL_alphaxgamma*(L'*b);
    
end

end

