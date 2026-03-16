function [CD, CDn] = getColumnDependency(A,B)
% returns the (normalized) column dependency (=the sum of scalar products
% of all normalized column vectors for rows i != j) for matrix A (with B)

% normalize columns
An = normcols(A);
% get size of matrix
[~,c] = size(A);
if nargin == 1
    % calculate column dependency
    CD = sum(sum(abs(triu(An'*An)))) - c;
    % calculate normalized column dependency
    CDn = 2*CD/(c*(c-1));
else
    % normalize columns
    Bn = normcols(B);
    % calculate column dependency
    CD = sum(sum(abs(triu(An'*Bn))));
    % calculate normalized column dependency
    CDn = 2*CD/(c*(c+1));
end
end

