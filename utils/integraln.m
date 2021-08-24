function int = integraln(fun, lb, ub, max_quad_points)

    assert(size(lb, 2) == 1)
    assert(size(ub, 2) == 1)
    assert(isa(fun, 'function_handle'))
    m = length(lb);
    if m >= 4
        assert(all(lb > -inf), 'If ndim>3 simple trapezoidal integration is used, which does not support inf bounds.')
        assert(all(ub < inf), 'If ndim>3 simple trapezoidal integration is used, which does not support inf bounds.')
    end
    xtest = fun([lb, ub]);
    assert(all(size(xtest) == [1, 2]))
    
    if nargin < 4
        max_quad_points = 1e5;
    end

    if m == 1
        int = integral(fun, lb, ub);
    elseif m == 2
        if any(isinf([lb(:); ub(:)]))
            int = integral2(@(X, Y) reshape(fun([X(:)'; Y(:)']), size(X)), lb(1), ub(1), lb(2), ub(2));
        else
            int = quad2d(@(X, Y) reshape(fun([X(:)'; Y(:)']), size(X)), lb(1), ub(1), lb(2), ub(2), ...
                'MaxFunEvals', max_quad_points);
        end
    elseif m == 3
        int = integral3(fun, lb(1), ub(1), lb(2), ub(2), pdf_lb(3), ub(3));
    else
        % Resort to simple trapezoidal integration.
        % How fine should the grid be?
        num_quad_points_per_dim = floor(max_quad_points^(1/m)-1);

        % create grid for integration
        trapz_grid_vecs = cell(1, m);
        trapz_grid_step = (ub - lb) ./ (num_quad_points_per_dim - 1);
        for ii = 1:m
            trapz_grid_vecs{ii} = lb(ii) + (0:num_quad_points_per_dim) * trapz_grid_step(ii);
        end
        trapz_grid_matrices = cell(1, m);
        [trapz_grid_matrices{:}] = ndgrid(trapz_grid_vecs{:});

        % evaluate pdf on grid
        trapz_grid_matrices_unrolled = ...
            cellfun(@(X) reshape(X, 1, []), trapz_grid_matrices, 'UniformOutput', false);
        Xq = cat(1, trapz_grid_matrices_unrolled{:});
        fun_vals = fun(Xq);

        % we need the results in the shape of an n-d array
        fun_vals_array_shape = cellfun(@length, trapz_grid_vecs, 'UniformOutput', false);
        fun_vals_array = reshape(fun_vals, fun_vals_array_shape{:});

        % perform trapezoidal integration
        int = trapezoidal_rule_nd_integral(trapz_grid_vecs, fun_vals_array, m);
    end
end