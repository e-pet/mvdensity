function eps = find_rbf_eps(rbf_centers, method, varargin)

    %% Input handling
    p = inputParser;

    % Positional arguments
    addRequired(p, 'rbf_centers', @ismatrix);
    addRequired(p, 'method', @isvector);
    addOptional(p, 'debug_info', false, @islogical);
    % Name-value arguments
    addParameter(p, 'rbf', [], @isfunction);
    addParameter(p, 'nodes', [], @ismatrix);
    addParameter(p, 'nodevals', [], @ismatrix);
    addParameter(p, 'center_eps_scale', [], @ismatrix);

    parse(p, rbf_centers, method, varargin{:})
    struct2vars(p.Results);

    %% Actual eps determination
    
    % Heuristics for epsmin and epsmax are mine
    center_DM = squareform(pdist(rbf_centers));
    epsmin = 1/max(center_DM(:));
    epsmax = 10/mean(center_DM(:));    
    
    if debug_info
        fprintf('Choosing RBF eps using %s.\n', method);
        disp(['Epsmin ', num2str(epsmin), ', Epsmax ', num2str(epsmax)]);    
    end
    
    if strcmp(method, 'loocv') || strcmp(method, 'auto')
        % Find optimal shape parameter epsilon using leave-one-out cross-
        % validation (LOOCV) as proposed by Rippa (1999):
        % "An algorithm for selecting a good value for the parameter c in radial basis
        % function interpolation", and as described in detail in Fasshauer/Zhang
        % (2007): "On Choosing “Optimal” Shape Parameters for RBF Approximation"
        % https://link.springer.com/content/pdf/10.1007/s11075-007-9072-8.pdf
        DM = pdist2(nodes, rbf_centers);
        [eps, ~] = fminbnd(@(eps) CostEps(eps, rbf, DM, nodevals), epsmin, epsmax);
        
        if strcmp(method, 'auto') && (abs(epsmax-eps)/epsmax < 1e-3 || abs(eps-epsmin)/epsmin < 1e-3)
            % For some reason, the optimization apparently does not work well.
            % Use Hardy's rule as fall-back.
            if debug_info
                fprintf('LOOCV eps %f is out of bounds. Using Hardy''s heuristic formula as fall-back.', eps);
            end
            eps = eps_hardy(center_DM);
        end
        
    elseif strcmp(method, 'hardy')
        eps = eps_hardy(center_DM);
        
    elseif strcmp(method, 'random')
        % Choose eps randomly within a reasonable range as some people have
        % shown that using different eps values can increase numerical
        % stability.
        eps = epsmin + rand(1, size(center_DM, 1)) * (epsmax - epsmin);
        
    elseif strcmp(method, 'prop')
        % Choose varying eps proportional to the user-provided vector
        % node-eps-scale, between epsmin and epsmax.
        eps_scale = (center_eps_scale - min(center_eps_scale)) / (max(center_eps_scale) - min(center_eps_scale));
        assert(max(eps_scale) == 1);
        assert(min(eps_scale) == 0);
        eps = epsmin + eps_scale * (epsmax - epsmin);
    else
        error('Unknown method requested.');
    end
    
    if debug_info
        fprintf(['Final selected eps: ', num2str(eps(1:min(length(eps(:)), 10))')]);
    end
end


% Define LOOCV cost function to minimize; code (almost completely) copy-pasted from Fasshauer/Zhang (2007)
function ceps = CostEps(eps, rbf, DM, rhs)
    A = rbf(DM, eps);
    invA = pinv(A);
    errorvector=(invA * rhs) ./ diag(invA);
    ceps = norm(errorvector);
end


function eps = eps_hardy(center_DM)
    % Use the rule of thumb proposed by Hardy, R.L.: Multiquadric equations 
    % of topography and other irregular surfaces. J. Geophys. Res. 76, 
    % pp. 1905–1915 (1971).
    % Construct distance matrix without 0 diagonal
    N = size(center_DM, 1);
    DMd = center_DM;
    DMd(DMd == 0) =[];
    DMd = reshape(DMd, N, N-1);
    % average of distance to closest neighbor
    d = mean(min(DMd, [], 2));
    eps = 1/(0.815*d);
end