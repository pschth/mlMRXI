function kappa_fro = frocond( A )
% returns the Frobenius condition number of matrix A

kappa_fro = norm(A, 'fro') * norm( pinv(A), 'fro');

end

