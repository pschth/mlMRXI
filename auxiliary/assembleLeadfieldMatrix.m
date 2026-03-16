function L_new = assembleLeadfieldMatrix( L_dict, I )
% assembles leadfield matrix L_new from activation current matrix I using
% the dictionary leadfield matrix L_dict with nCoils individual coil
% activations with unit current (i.e. I_dict = eye(nCoils))

% INPUT
% L_dict ... nCoils*nSensors x nVox; Dictionary leadfield matrix containing
% the leadfield matrices of all coils in use with individual activation and
% concatenated below each other
% I ... nCoils x nActivations; contains the coil currents

% OUTPUT
% L_new ... nActivations*nSensors x nVox; Leadfield matrix according to
% activation currents I

nCoils = size(I,1);
nSensors = size(L_dict,1)/nCoils;

if nSensors/round(nSensors) ~= 1
    error('nSensors no integer: Number of rows of L_dict divided by number of coils (number of rows of I) is no integer.');
end

Lv = vectorizeLeadfieldMatrix(L_dict, 'f', nSensors);
L_new = vectorizeLeadfieldMatrix(Lv*I, 'b', nSensors);

end

