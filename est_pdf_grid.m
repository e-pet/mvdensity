function pdf_handle = est_pdf_grid(data, nbin, exclusion_quantile, smoothing_method, debug_info, data_labels)
% EST_PDF_GRID Estimate a multivariate PDF by histogram interpolation
%
% This method simply covers m-dimensional space with a grid of nbins along each dimension, counts samples
% in each bin, and then performs either smoothing (if smoothing_method='rbf') or interpolation (else) of the
% histogram count values. To prevent negative values, the estimation is performed on a nonlinear transformation
% of the count values, which ensures positivity. The method returns a function handle to the estimated PDF 
% surface, which is vectorized and can be efficiently evaluated at many datapoints.
%
% data: samples from which to identify the PDF. Can be either Nxm or mxN, where N=number of samples and
%          m=number of dimensions.
% nbin [optional]: number of bins to use. The same number of bins is used along each dimension, so in total
%                        there will be m^nbin regions. If not provided, nbin is determined automatically using a 
%                        heuristic formula.
% exclusion_quantile [optional]: Data samples below the q-quantile and above the (1-q)-quantile, where
%                                            q=exclusion_quantile, will be fully excluded from the analysis. This can be
%                                            used in order to prevent outliers from distorting histogram bin placement.
% smoothing_method [optional]: If 'rbf', radial basis function (RBF)-based smoothing is used. If anything else,
%                                             the corresponding "griddedInterpolant" method is used. The default choice is
%                                             makima interpolation.
% debug_info [optional]: Should debug plots be generated and debug info be printed to the command line?
%                                 Default: false.
% data_labels [optional]: If provided, must be a cell array of length m containing variable labels, which will be
%                                  used to label some of the plots generated if debug_info==true.
%
% RETURNS: pdf_handle, a function handle that accepts a single argument of dimension Nqxm, where Nq is the
%                 number of query points at which the estimated PDF is to be evaluated.
%
% Eike Petersen, August 2021.
%

    [N, m] = size(data);
    if m > N
        data = data';
        [N, m] = size(data);
    end
    
    if nargin < 2 || isempty(nbin)
        % Heuristic for choosing the number of regions. Note that we choose the number of regions / 
        % tiles equal in each dimension, because I am not aware of any clever way to do this more
        % efficiently.
        % Rationale for the formula: in 1-D, I consider ~10 bins to be reasonable @ 1000 samples,
        % and ~50 regions @ 100,000 samples. Also, nbin should grow _significantly sublinearly
        % with nsample, so we choose nbin ~ log(nsample).
        % In particular, in 1-D, we choose nbin = a + b log(nsample), and by inserting the two
        % above datapoints, we obtain a=-50 and b=8.6859.
        % Now, how to extend this to m-D? Of course, nregion = nbin^m. So what I first did is
        % nbin = (a + b log(nsample))^(1/m). But that grows _really_ slowly in higher dimensions,
        % e.g., in 5-D we still only have three histogram bins @ 1e6 datapoints. So I chose to scale
        % it with sqrt(m) instead of m, yielding the below formula.
        % Finally, the number of bins should always be at least 3, and the number of regions should
        % always be below 1 million.
        nbin = min(max(3, floor(max((8.6859*log(N)-50), 0)^(1/sqrt(m)))), floor(1e6^(1/m)));
    end
    
    if nargin < 3 || isempty(exclusion_quantile)
        exclusion_quantile = min(0.01, 10/N);
    end
    
    if nargin < 4 || isempty(smoothing_method)
        smoothing_method = 'makima';
    end    
    
    if nargin < 5
        debug_info = false;
    end
    
    if nargin < 6
        data_labels = [];
    end
    
    data_range = max(data) - min(data);
    
    % Scale to standard deviation 1.
    data_scale = std(data, 'omitnan');
    data = data ./ data_scale;

    
    %% A. Construct an n-dimensional histogram
    
    % Select the range to cover with tiles in each dimension. Everything
    % outside these bounds will be ignored completely, because taking it
    % into account in the outer tile would lead to some distortion of the 
    % density in that tile.
    % Note that we do not want to simply select the min-max range, because
    % in the presence of strong outlier, this would dilute the regions
    % strongly over the data space.
    sd_minmax = nanquantile(data, [exclusion_quantile, 1-exclusion_quantile]);
    sd_min = sd_minmax(1, :);
    sd_max = sd_minmax(2, :);
    sd_range = sd_max - sd_min;
    sd_step = sd_range / nbin;
    
    grid_lb = sd_min .* data_scale;
    grid_ub = sd_max .* data_scale;
    
    pdf_lb = (sd_min-0.5*sd_step) .* data_scale;
    pdf_ub = (sd_max+0.5*sd_step) .* data_scale;
    is_within_bounds = @(X) all((X > pdf_lb) & (X < pdf_ub), 2);
    
    % drop everything outside the selected data range
    data_reduced = data(all(data > sd_min & data < sd_max, 2), :);
    
    % Build a matrix containing all regular grid positions.
    grid_vecs = cell(1, m);
    grid_vec_idces = cell(1, m);
    bin_bounds = cell(1, m);
    for ii = 1:m
        grid_vecs{ii} = sd_min(ii) + (-0.5:1:nbin+0.5)*sd_step(ii);
        grid_vec_idces{ii} = 1:nbin+2;
        bin_bounds{ii} = [sd_min(ii) + (0:1:nbin)*sd_step(ii)];
    end
    grid_matrices = cell(1, m);
    [grid_matrices{:}] = ndgrid(grid_vecs{:});
    grid_idces_matrices = cell(1, m);
    [grid_idces_matrices{:}] = ndgrid(grid_vec_idces{:});
    grid_idces = zeros((nbin+2)^m, m);
    for ii = 1:m
        grid_idces(:, ii) = grid_idces_matrices{ii}(:);
    end

    [bin_flags, bin_counts, bin_flag_idx] = hist_bin_counts(data_reduced, bin_bounds{:});
    [mf, ~] = size(bin_flags);
    
    all_flag_combinations = unique(nchoosek(repmat(1:nbin+2, 1, m), m), 'rows');
    all_flag_counts = 1e-5 * min(bin_counts) * ones(size(all_flag_combinations, 1), 1);
    all_tile_centers = sd_min - 1.5 * sd_step + all_flag_combinations .* sd_step;
    
    for ii = 1:mf
        combination_mask = all(all_flag_combinations == bin_flags(ii, :) + 1, 2);
        assert(sum(combination_mask) == 1);
        all_flag_counts(combination_mask) = bin_counts(ii);
        all_tile_centers(combination_mask, :) = mean(data_reduced(bin_flag_idx == ii, :));
    end
    assert(sum(all_flag_counts >= min(bin_counts)) == mf);
    
    % Set up m-dimensional array of count values
    if m >= 2
        count_array_shape = cellfun(@length, grid_vecs, 'UniformOutput', false);
        count_array = zeros(cell2mat(count_array_shape));
        for ii = 1:length(all_flag_counts)
            tile_ii_idces = num2cell(all_flag_combinations(ii, :));
            count_array(tile_ii_idces{:}) = all_flag_counts(ii);
        end
    else
        count_array = zeros(size(all_flag_counts));
        for ii = 1:length(count_array)
            count_array(all_flag_combinations(ii)) = all_flag_counts(ii);
        end
    end
    
    %% B. Interpolate histogram counts
    % Now obtain approximate "counts" for each data sample by interpolating
    % the n-dimensional histogram.
    % To guarantee positivity (and maybe even make the interpolation task
    % easier?), do this in the softplus^-1 space: transform the data by
    % log(exp(y)-1) before interpolation. This justifies transforming the
    % resulting values by log(1+exp(x)), the so-called soft-plus function,
    % a smooth approximation to the "rectifier" max(0, x), which guarantees
    % positivity of the resulting values. 
    % The soft-plus function has the added benefit that it behaves linearly
    % for large values and thusly reacts less strong to outliers compared
    % to the log/exp transformation pair (which would also guarantee
    % positivity).

    if strcmp(smoothing_method, 'rbf')
        % Using RBFs, we can use irregularly scattered data. Therefore, we
        % use the tile centers of mass as the interpolation points.
        % Since we want to smooth - not interpolate - the histogram, select
        % significant tiles to use. Significance is defined here as
        % abs(L/std(L)), where L denotes the (discrete) Laplacian of
        % the data counts. In addition to the tiles found to be significant,
        % we also add all corner tiles to ensure nice behavior towards the
        % edges.

        % Approximate Laplacian over the histogram. Assume equal spacing (this
        % is equivalent (??) to a rescaling).
        L = del2(count_array);

        % Calculate bin significance following Eq. (16) in Allison (1993):
        % Multiquadric Radial Basis Functions for Representing Multidimensional
        % High Energy Physics Data.
        S = abs(L ./ sqrt(var_of_Laplacian(count_array)));        
        
        if nbin^m <= 5
            target_num_sig_tiles = nbin^m;
        else
            target_num_sig_tiles = 5 + (nbin^m - 5)^(log(30)/log(75));
        end
        sig_cutoff = quantile(S(:), 1-target_num_sig_tiles/(nbin^m));
        significant_arr_mask = S >= sig_cutoff;

        sig_idces = cell(m, 1);
        [sig_idces{:}] = ind2sub(size(significant_arr_mask), find(significant_arr_mask(:)));
        significant_grid_idces = horzcat(sig_idces{:});
        n_sig_tiles = sum(significant_arr_mask, 'all');
        assert(all(size(significant_grid_idces) == [n_sig_tiles, m]));
        
        boundary_flag = all_tile_centers < sd_min | all_tile_centers > sd_max;
        corner_rows = all(boundary_flag, 2);
        n_corner_tiles = sum(corner_rows);
        
        rbf_centers = zeros(n_sig_tiles + n_corner_tiles + 1, m);
        for ii = 1:n_sig_tiles
            rbf_centers(ii, :) = all_tile_centers(all(significant_grid_idces(ii, :) == all_flag_combinations, 2), :);
        end
        % Significance of each center
        rbf_center_sig = zeros(n_sig_tiles + n_corner_tiles + 1, 1);
        rbf_center_sig(1:n_sig_tiles) = S(sub2ind(size(S), sig_idces{:}));
        assert(all(rbf_center_sig(1:n_sig_tiles) >= sig_cutoff));
        % Add all corners
        rbf_centers(n_sig_tiles+1:end-1, :) = all_tile_centers(corner_rows, :);
        corner_idces = mat2cell(grid_idces(corner_rows, :), sum(corner_rows), ones(m, 1));
        rbf_center_sig(n_sig_tiles+1:end-1) = S(sub2ind(size(S), corner_idces{:}));
        % Add maximum peak
        [~, peakidx] = max(all_flag_counts);
        peak_nd_idx = mat2cell(grid_idces(peakidx, :), 1, ones(m, 1));
        rbf_centers(end, :) = all_tile_centers(peakidx, :);
        rbf_center_sig(end, :) = S(sub2ind(size(S), peak_nd_idx{:}));
        % Remove duplicates (a corner or the peak may also be considered significant)
        [rbf_centers, unique_idces, ~] = unique(rbf_centers, 'rows');
        rbf_center_sig = rbf_center_sig(unique_idces);
        
        % Choose variable RBF eps proportional to bin significance
        eps = find_rbf_eps(rbf_centers, 'prop', true, 'center_eps_scale', -rbf_center_sig);
        interp_func = interp_rbf(all_tile_centers, soft_plus_inv(all_flag_counts), rbf_centers, 'multiquadric', eps, ...
            debug_info);  %, min(1e-2, 10^(-log10(N))));
        pdf_handle_unnormalized = @(Xq) is_within_bounds(Xq) .* max(soft_plus(interp_func(Xq ./ data_scale)), 0);
    else    
        % Interpolation only works with regularly gridded data, whence we use the regular histogram grid points as 
		% interpolation points.
        % See https://blogs.mathworks.com/cleve/2019/04/29/makima-piecewise-cubic-interpolation/
        % for a nice explanation and some examples of what makima does, which is the default choice.
        interpfun = griddedInterpolant(grid_vecs, soft_plus_inv(count_array), smoothing_method, 'nearest'); % 'nearest' is ExtrapolationMethod
        pdf_handle_unnormalized = @(Xq) is_within_bounds(Xq) .* soft_plus(interpfun(Xq ./ data_scale));
    end
    
    % Calculate an approximate integral over the estimated unnormalized pdf
    % to provide approximate normalization of the resulting PDF estimate.
    % This is only an approximation of the actual normalizing constant, which
    % is all the sadder as that constant could be calculated exactly relatively
    % easily: one simply needs to calculate the integral over all the RBFs / all
    % the cubic splines used for the interpolation. This is somewhat
    % tedious in N dimensions, though...
    warning('off', 'MATLAB:quad2d:maxFunEvalsFail');
    warning('off', 'MATLAB:quad2d:maxFunEvalsPass');
    approx_integral = integraln(@(x) pdf_handle_unnormalized(x')', pdf_lb', pdf_ub');
    pdf_handle = @(Xq) pdf_handle_unnormalized(Xq) ./ approx_integral;
    
    
    %% C. Debug ouptut and plots
    
    if debug_info
        disp('---------------------------');
        disp(['Histogram-based PDF estimation using ', smoothing_method]);
        disp(['Num Data Points: ', num2str(N)]);
        disp(['Num Dimensions: ', num2str(m)]);
        disp(['SD Range: ', num2str(sd_range)]);
        disp(['Num Regions: ', num2str(nbin)]);
        disp(['Num Tiles: ', num2str((nbin+2)^m)]);
        disp(['Lower bound of region spanned by (real, non-artificial) histogram bins:\n  ', num2str(grid_lb)]);
        disp(['Upper bound of region spanned by (real, non-artificial) histogram bins:\n  ', num2str(grid_ub)]);        
        if strcmp(smoothing_method, 'rbf')
            disp(['Num RBF Centers: ', num2str(size(rbf_centers, 1))]);
        end
        
        fig_title = sprintf('PDF hist interpolation using %s', smoothing_method);
        fig1 = figure('NumberTitle', 'off', 'Name', fig_title);
        [~, AX, BigAx, ~, ~] = plotmatrix(data .* data_scale);
        if ~isempty(data_labels)
            for jj = 1:m
                ylabel(AX(jj, 1), data_labels{jj});
                xlabel(AX(m, jj), data_labels{jj});
            end
        end

        for ii = 1:m
            hold(AX(ii, ii), 'on');
            % Add tile lines to histogram of data dimension ii
            for jj = 0:nbin
                xline(AX(ii, ii), (sd_min(ii) + jj*sd_step(ii)) * data_scale(ii), '--');
            end
            
            % Add tile lines to all scatter plots in this plot row
            for jj = 1:m
                if ii~=jj
                    hold(AX(ii, jj), 'on');
                    for kk = 0:nbin
                        yline(AX(ii, jj), (sd_min(ii) + kk*sd_step(ii)) * data_scale(ii), '--');
                        xline(AX(ii, jj), (sd_min(jj) + kk*sd_step(jj)) * data_scale(jj), '--');
                    end
                end
            end
        end
        title(BigAx, fig_title);
        
        if m < 3
            fig2 = figure('NumberTitle', 'off', 'Name', fig_title);
            
            if strcmp(smoothing_method, 'rbf')
                h1 = subplot(1, 3, 1);
            else
                h1 = subplot(1, 2, 1);
            end
            
            if m == 1
                histogram(data .* data_scale);
                xlabel('x');
                ylabel('counts');
            elseif m == 2
                x1q = (sd_min(1):.25:sd_min(1)+sd_range(1)) * data_scale(1);
                x2q = (sd_min(2):.25:sd_min(2)+sd_range(2)) * data_scale(2);
                histogram2(data(:, 1) * data_scale(1), data(:, 2) * data_scale(2));
                xlim([min(x1q), max(x1q)]);
                ylim([min(x2q), max(x2q)]);
                xlabel('x1');
                ylabel('x2');
                zlabel('counts');
            end
            
            if strcmp(smoothing_method, 'rbf')
                h2 = subplot(1, 3, 2);
            else
                h2 = subplot(1, 2, 2);
            end
            
            if m == 1
                scatter(all_tile_centers * data_scale, all_flag_counts, 'filled');
                hold on;
                for ii = 0:nbin
                    xline((sd_min + ii*sd_step)*data_scale, '--');
                end
                xq = linspace(min(data .* data_scale) - 0.05*data_range, max(data .* data_scale) + 0.05*data_range, 100);
                vq = pdf_handle(xq');
                plot(xq, approx_integral .* vq);
                linkaxes([h1, h2], 'x');
                xlabel('x');
                ylabel('pdf');
            elseif m == 2
                [X1, X2] = meshgrid(x1q, x2q);
                vq = pdf_handle([X1(:), X2(:)]); 
                surf(X1, X2, approx_integral .* reshape(vq, size(X1)), 'FaceAlpha', '0.5', 'EdgeColor', 'black', 'FaceColor', '#D95319');
                hold on;
                scatter3(all_tile_centers(:, 1) * data_scale(1), all_tile_centers(:, 2) * data_scale(2), ...
                    all_flag_counts, 'filled');
                xlabel('x1');
                ylabel('x2');
                zlabel('pdf');
            end

            if strcmp(smoothing_method, 'rbf')
                h3 = subplot(1, 3, 3);
                if m == 1
                    vqS = interp1(grid_vecs{1}, S, xq/data_scale, 'nearest');
                    plot(xq, vqS);
                    for ii = 0:nbin
                        xline((sd_min + ii*sd_step)*data_scale, '--');
                    end                
                    linkaxes([h1, h2, h3], 'x');
                    xlim([min(xq), max(xq)]);
                    xlabel('x');
                    ylabel('Bin Significance = L/std(L)');
                elseif m == 2
                    Vq = interp2(grid_matrices{1}', grid_matrices{2}', S', X1/data_scale(1), X2/data_scale(2), 'nearest'); 
                    surf(X1, X2, Vq, 'FaceAlpha', '0.5', 'EdgeColor', 'black', 'FaceColor', '#D95319');
                    xlabel('x1');
                    ylabel('x2');
                    zlabel('Bin Significance = L/std(L)');
                end
            end
            gridTitle(fig2, fig_title);
        end
    end
end


function Laplacian_filter_arr = discrete_Laplacian(m)
    % Set up m-dimensional array with size 3 in each dimension
    Laplacian_filter_arr = zeros([repmat(3, 1, m), 1]);
    % Set all entries where some dimension has index 2 to 1
    dim_vecs = mat2cell(repmat(1:3, m, 1), ones(m, 1));
    idces = cell(m, 1);
    [idces{:}] = ndgrid(dim_vecs{:});
    idces_arr = reshape(cat(m+1, idces{:}), [], m);
    idx_with_one_2_flag = sum(idces_arr == 2, 2) == m-1;
    dim_positions = mat2cell(idces_arr(idx_with_one_2_flag, :), sum(idx_with_one_2_flag), ones(1, m));
    lin_array_idx_with_2 = sub2ind(size(Laplacian_filter_arr), dim_positions{:});
    Laplacian_filter_arr(lin_array_idx_with_2) = 1;
    % Set the center element to -2*m
    center_pos = mat2cell(repmat(2, m, 1), ones(m, 1));
    Laplacian_filter_arr(center_pos{:}) = -2*m;
end

function VL = var_of_Laplacian(arr)
    m = ndims(arr);
    var_L = discrete_Laplacian(m).^2;
    VL = convn(arr, var_L, 'same');
end