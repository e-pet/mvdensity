function interp = interp_rbf(nodes, nodevals, centers, kernel, eps, debug_info, ridge_k)
% Perform N-D interpolation and smoothing using radial basis functions (RBFs)
%
% Efficient procedure to perform n-dimensional interpolation and smoothing. 
% Data can be scattered arbitrarily and do not need to conform to a regular 
% grid. Can also be used for extrapolation (with reduced precision, obviously).
% 
% See, e.g., https://en.wikipedia.org/wiki/Radial_basis_function_interpolation
%
% INPUTS:
%   - "nodes" should be an Nxp matrix, where N is the number of nodes and p
%   is the number of dimensions, i.e., each row in this matrix specifies
%   the coordinates of one measured node.
%   - "nodevals" should be a vector of length N which contains the desired
%   function values at each node.
%   - "centers" should be a vector of length M<=N which contains the
%   locations of the RBF center locations. These may be identical to the
%   nodes, but they do not have to. M==N results in interpolation being
%   performed, M<N results in a smoothed estimate.
%   - "kernel" (optional) is a character vector specifying the RBF kernel to  
%   be used. Currently, only an "inverse multiquadric" basis function is 
%   implemented, but others could be added easily.
%   - "eps" (optional) specifies the value of the shape parameter in the
%   kernel function (if there is one). If not provided, eps is selected
%   automatically by means of leave-one-out cross-validation (LOOCV).
%   - "debug_info" (optional) specifies whether detailed information should
%   be displayed for debugging purposes.
%   - "ridge_k" (optional) can be used to specify the ridge regression
%   regularization constant. Default is 1e-5. Check the Matlab "ridge"
%   documentation to see the effect of this constant.
%
% RETURNS: A function handle that takes query nodes and returns
% interpolated function values at those nodes. "qnodes" should be an Mxp 
% matrix, where M is the number of query nodes and p is the number of 
% dimensions, i.e., each row in this matrix specifies the coordinates of
% one query node.
%
% COMMENT: RBF interpolation works nicely and easily BUT cannot easily
% guarantee positivity of the resulting interpolation. Even when both
% data and RBF functions are positive everywhere, weights may be
% negative to ensure meeting some data points exactly. Preventing this
% seems to be a non-trivial endavour, see, e.g., 
% https://aip.scitation.org/doi/abs/10.1063/1.4980915?journalCode=apc
% A simple hack to achieve positivity is the take the log of the
% weights before interpolation, and then transforming back after
% transformation.
%
% EXAMPLE: Recreate sin*cos from 25 randomly sampled data points
%
% nodes = rand(25, 2);
% nodevals = sin(2*pi*nodes(:, 1)) .* cos(2*pi*nodes(:, 2));
% [X, Y] = meshgrid(0:.02:1, 0:.02:1);
% interp_func = interp_rbf(nodes, nodevals);
% gridvals = interp_func([X(:), Y(:)]);
% figure;
% scatter3(nodes(:, 1), nodes(:, 2), nodevals);
% hold on;
% surf(X, Y, reshape(gridvals, size(X)), 'FaceAlpha', 0.5, 'EdgeColor', 'none')
%

[N, mn] = size(nodes);
[M, mc] = size(centers);
assert(length(nodevals) == N);
%assert(M <= N);
assert(mn == mc);

% Calculate distance matrix
DM = pdist2(nodes, centers);

% Define RBF kernel function
if nargin < 4 || isempty(kernel) || strcmp(kernel, 'invmultiquadric')
    rbf = @(dist, eps) 1./sqrt(1 + repmat(eps', size(dist, 1), 1).^2 .* dist.^2);
elseif strcmp(kernel, 'multiquadric')
    rbf = @(dist, eps) sqrt(1 + repmat(eps', size(dist, 1), 1).^2 .* dist.^2);
else
    error('Requested kernel function is unknown or not implemented.')
end

if nargin < 5 || isempty(eps) || strcmp(eps, 'auto')
    eps = find_rbf_eps(nodes, nodevals(:), centers, rbf, debug_info);
end

assert(size(eps, 2) == 1)
if size(eps, 1) == 1
    eps = repmat(eps, M, 1);
end

if nargin < 6 || isempty(debug_info)
    debug_info = false;
end

if nargin < 7
    ridge_k = 1e-5;
end

% Calculate Nx1 RBF weights by solving a linear NxN system
A = rbf(DM, eps);
if debug_info
    disp(['RBF cond: ', num2str(cond(A))]);
end
% opts.SYM = true; opts.POSDEF = false;
%rbf_weights = linsolve(A, nodevals(:), opts);
rbf_weights = ridge(nodevals(:), A, ridge_k, 0);
assert(all(size(rbf_weights) == [M+1, 1]));

interp = @(qnodes) rbf_interpolant(centers, rbf_weights, rbf, eps, qnodes);

end


function qvals = rbf_interpolant(centers, rbf_weights, rbf, eps, qnodes)
    [N, p] = size(centers);
    assert(size(qnodes, 2) == p);
    [M, ~] = size(qnodes);
    % Construct MxN distance matrix between RBF nodes and query nodes. Element
    % (i,j) contains the distance between query node i and RBF node j.
    warning('off', 'all'); % there is a double-to-single conversion warning in here when this is called from fsurf
    DMq = pdist2(qnodes, centers);
    warning('on', 'all')
    assert(all(size(DMq) == [M, N]));

    % Calculate interpolated function values at query points
    qvals = rbf_weights(1) + rbf(DMq, eps) * rbf_weights(2:end);
    assert(all(size(qvals) == [M, 1]));
end