function  A_normeddim  = normdim( A, dim )
% normalizes provided dimension dim of matrix A

if nargin < 2
    dim = 1;
elseif isempty(dim)
    dim = 1;
end

A_normeddim = bsxfun(@rdivide, A, rssq(A,dim));
if any(isnan(A_normeddim(:)))
    warning('Dimensional normalization resulted in NaN values.');
end
end
