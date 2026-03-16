function  A_colnormed  = normcols( A )
% normalizes the columns (dimension: 1) of matrix A

% divide columns of A by norm of columns of A
A_colnormed = bsxfun(@rdivide, A, rssq(A,1));

% return a warning if NaN values are detected (most likely because of
% division by 0)
if any(isnan(A_colnormed(:)))
    warning('Column normalization resulted in NaN values.');
end
end
