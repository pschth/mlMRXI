function Aout = vectorizeLeadfieldMatrix( A, FoB, nSens )
% FoB = 'f': Rearranges the nSensors x nVoxels blocks of a leadfield matrix
% into nActivation column matrices.
% FoB = 'b': Reverts the rearranged, vectorized leadfield matrix to normal.
% 
% INPUT
% A - FoB: 
% 'f' nSensors*nActivations x nVoxels matrix; leadfield matrix
% 'b' nSensors*nVoxels x nActivations matrix; vectorized leadfield matrix

% FoB - string; 
% 'f': reorganizes each leadfield matrix of a single activation into the column
% vectors of a vectorized leadfield matrix Aout.
% 'b': reorganizes a vectorized leadfield matrix A back into the
% original leadfield matrix Aout.
% 
% nSens - scalar; number of sensors
% 
% OUTPUT
% Aout - FoB: 
% 'f' nSensors*nVoxels x nActivations matrix; vectorized leadfield matrix
% 'b' nSensors*nActivations x nVoxels matrix; leadfield matrix


if FoB == 'f'
    nVox = size(A,2);
    nAct = size(A,1)/nSens;

    if ~(nSens==round(nSens))
        error('Number of leadfield matrix rows is not equal to nActivations*nSensors.');
    end

    Aout = reshape( permute( reshape( A', [nVox nSens nAct]), [2 1 3] ), [nVox*nSens nAct]);

elseif FoB == 'b'
    nAct = size(A,2);
    nVox = size(A,1)/nSens;
    
    if ~(nVox==round(nVox))
        error('Number of vectorized leadfield matrix rows is not equal to nVoxels*nSensors.');
    end
    
    Aout = reshape( permute( reshape(A, [nSens nVox nAct]), [1 3 2]), [nSens*nAct nVox]);
    
else
    error('FoB can only be "f" or "b".');
end

end

