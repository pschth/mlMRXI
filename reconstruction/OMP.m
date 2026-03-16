function x = OMP(A, b, m)
% performs orthogonal matching pursuit
% INPUT:
% A ... dictionary matrix
% b ... measurement vector
% m ... desired sparsity level
% OUTPUT:
% x ... sparse data vector (A*x = b) with sparsity level m

% initialize indices variable
T = false(size(A,2),1);

% initialize residual
b = b(:);
r = b;
for ii = 1:m
    % find maximum element of A'*b
    [~,i] = max(abs(A'*r));
    % store index
    T(i(1)) = true;
    % solve least squares problem to get estimate of atom amplitude
    xx = A(:, T)\b;
    % calculate new approximation of the data and the new residual
    a = A(:, T) * xx;
    r = (b - a);
end
x = zeros(size(A,2),1);
x(T) = xx;
end

