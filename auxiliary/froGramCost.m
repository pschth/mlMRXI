function cost = froGramCost( L, I, nSensors )
% returns the Frobenius norm of the Gram matrix of the leadfield matrix
% calculated from the vectorized leadfield matrix Lv times the coil current
% matrix I.
% 
% INPUT
% L - nSensors*nActivations x nVoxel OR nSensors*nVoxels x nCoils matrix
% (vectorized); standard OR vectorized (dictionary) leadfield matrix (see
% function vectorizeLeadfieldMatrix.m for description)  
% 
% I (if L is vectorized) - nCoils x nActivations matrix; coil current matrix
% 
% nSensors (if L is vectorized) - scalar; number of sensors
% 
% OUTPUT
% cost - scalar; see description of function

if nargin > 1
    if size(L,2) ~= size(I,1)
        error('Matrix sizes of Lv and/or I do not fit.');
    end
    % calculate actual leadfield matrix
    L = vectorizeLeadfieldMatrix(L*I,'b',nSensors);
end
% column-normalize leadfield matrix
L_curr_cn = normcols(L);
nVox = size(L,2);
% calculate Frobenius norm of Gram matrix
cost = norm(L_curr_cn'*L_curr_cn - eye(nVox),'fro');

end

