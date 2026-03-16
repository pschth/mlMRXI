function V = create3DVoxelDictionary( L, nSensors )
% creates n_Voxel dictionaries V with n_Sensors x n_Coils entries, so that
% the ith column of the resulting leadfield matrix is V_i*I with I being
% the coil current vector.
%
% Output: 
% V ... n_Sensors x n_Coils x n_Voxel matrix
%
% Input: 
% L ... n_Sensors*n_Coils x n_Voxel matrix; leadfield matrix
% nSensors ... scalar; number of sensors

if numel(nSensors) == 1
   nCoils = size(L,1)/nSensors;
   nVox = size(L,2);
   V = reshape(L,nSensors,nCoils,nVox);    
end


end

