function FP = getFramePotential(A,B)
% returns the squared row dependency (= frame potential; the sum of squared
% scalar products of all normalized row vectors for rows i != j) for matrix
% A (with B)

% normalize rows
A = normdim(A,2);
if nargin == 1
    % calculate frame potential
    FP = sum(sum(triu(A*A', 1).^2));
else
    % normalize rows
    B = normdim(B,2);
    % calculate frame potential
    FP = sum(sum((A*B').^2));
end
end


