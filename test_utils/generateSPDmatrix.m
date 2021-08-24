function A = generateSPDmatrix(n)
% Generate a dense n x n symmetric, positive definite matrix
% due to https://math.stackexchange.com/a/358092/171782
A = rand(n,n); % generate a random n x n matrix

% construct a symmetric matrix
A = 0.5*(A+A');

% since A(i,j) < 1 by construction and a symmetric diagonally dominant matrix
%   is symmetric positive definite, which can be ensured by adding nI
A = A + n*eye(n);

end