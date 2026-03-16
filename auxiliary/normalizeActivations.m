function [L_normed, froL_single] = normalizeActivations(L, nSensors)
% normalizes the individual activations of a system matrix L with respect
% to the Frobenius norm
nActs = size(L,1)/nSensors;
froL_single = zeros(nActs,1);
L_normed = zeros(size(L));
for a = 1:nActs
    si = (a-1)*nSensors+1:a*nSensors;
    froL_single(a) = norm(L(si,:), 'fro');
    L_normed(si,:) = L(si,:)./froL_single(a);
end

