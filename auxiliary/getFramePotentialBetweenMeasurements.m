function [FPM, FP, FP_matrix] = getFramePotentialBetweenMeasurements(A,nSensors)
% returns frame potential between all measurement matrices

% reshape A to nSensors x nVoxel x nActivations matrix and normalize rows
[nSA, nVoxel] = size(A);
nActs = nSA/nSensors;
A = permute(reshape(normcols(A'),nVoxel,nSensors,nActs),[2 1 3]);

FP_matrix = zeros(nActs,nActs);
for i = 1:nActs
    FP_matrix(i,i) = sum(sum(triu(A(:,:,i)*A(:,:,i)',1).^2));
    for j = i+1:nActs
        FP_matrix(i,j) = sum(sum((A(:,:,i)*A(:,:,j)').^2));
    end
end
FP = sum(sum(FP_matrix));
FPM = sum(sum(triu(FP_matrix,1)));
end


