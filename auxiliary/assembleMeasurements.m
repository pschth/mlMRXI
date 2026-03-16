function b_new = assembleMeasurements( b_dict, I )
% assembles measurements b_new from activation current matrix I using the
% dictionary measurement vector b_dict with nCoils individual coil
% activations with unit current (i.e. I_dict = eye(nCoils))

% INPUT
% b_dict ... nCoils*nSensors vector; Dictionary measurement vector
% containing the measurements of all coils in use with individual
% activation and concatenated below each other
% I ... nCoils x nActivations; contains the coil currents

% OUTPUT
% b_new ... nActivations*nSensors vector; measurement vector according to
% activation currents I

nCoils = size(I,1);
nSensors = length(b_dict)/nCoils;

if nSensors/round(nSensors) ~= 1
    error('nSensors no integer: Number of rows of b_dict divided by number of coils (number of rows of I) is no integer.');
end

b_mat = reshape(b_dict, nSensors, nCoils);
b_mat_new = b_mat*I;
b_new = b_mat_new(:);

end

