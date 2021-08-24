function performance = run_example(test_pdf, pdf_est_method, N, do_plot, plot_title, ...
                                   export_example_file_samples, export_example_file_pdf)

    xn = random(test_pdf, N);
    
    [~, m] = size(xn);
    estimated_pdf = pdf_est_method(xn, do_plot);

    if nargin > 5 && ~isempty(export_example_file_samples)
        writematrix([(1:N)' xn], export_example_file_samples);
    end    
    
    if nargin > 6 && ~isempty(export_example_file_pdf)
        if m == 1
            xq = linspace(min(xn), max(xn), 1000);
            truepdfvals = pdf(test_pdf, xq');
            estpdfvals = estimated_pdf(xq');
            writematrix([xq' truepdfvals estpdfvals], export_example_file_pdf);
        else
            error('PDF export is currently not implemented in multiple dimensions, but could be.');
        end
    end
    
    % Jensen-Shannon divergence
    warning('off', 'MATLAB:quad2d:maxFunEvalsFail');
    warning('off', 'MATLAB:quad2d:maxFunEvalsPass');
    performance = js_div_nd_cont(@(x) pdf(test_pdf, x), @(x) estimated_pdf(x), min(xn), max(xn));
    
    if m == 1 && do_plot
        fig = figure('NumberTitle', 'off', 'Name', sprintf('%s, N = %d', plot_title, N));
        h1 = subplot(1, 2, 1);
        histogram(xn);
        xlabel('x');
        ylabel('counts');
        h2 = subplot(1, 2, 2);
        fplot(@(x) pdf(test_pdf, x')', [min(xn), max(xn)]);
        hold on;
        fplot(@(x) estimated_pdf(x')', [min(xn), max(xn)]);
        leg = legend('Truth', 'Estimation', 'AutoUpdate', 'off');
        % Enable clicking on axis to hide lines
        leg.ItemHitFcn = @hileline;
        xlabel('x');
        ylabel('p(x)');
        linkaxes([h1, h2], 'x');
        gridTitle(fig, sprintf('N = %d', N));

    elseif m == 2 && do_plot
        fig = figure('NumberTitle', 'off', 'Name', sprintf('%s, N = %d', plot_title, N));
        subplot(1, 2, 1);
        histogram2(xn(:, 1), xn(:, 2));
        xlabel('x1');
        ylabel('x2');
        x1min = min(xn(:, 1));
        x1max = max(xn(:, 1));
        x2min = min(xn(:, 2));
        x2max = max(xn(:, 2));
        xlim([x1min, x1max]);
        ylim([x2min, x2max]);
        zlabel('counts');
        subplot(1, 2, 2);
        fsurf(@(X, Y) reshape(pdf(test_pdf, [X(:), Y(:)]), size(X)), [x1min x1max x2min x2max], ...
            'FaceAlpha', 0.5, 'EdgeColor', 'black', 'FaceColor', '#0072BD');
        hold on;
        warning('off', 'MATLAB:fplot:NotVectorized');  % TBH I don't know why this is thrown: the function IS vectorized
        fsurf(@(X, Y) reshape(estimated_pdf([X(:), Y(:)]), size(X)), [x1min x1max x2min x2max], ...
            'FaceAlpha', 0.5, 'EdgeColor', 'black', 'FaceColor', '#D95319');
        leg = legend('Truth', 'MAKIMA', 'Location', 'best', 'AutoUpdate', 'off');
        % Enable clicking on axis to hide lines
        leg.ItemHitFcn = @hileline;
        xlabel('x1');
        ylabel('x2');
        zlabel('p(x1,x2)');
        gridTitle(fig, sprintf('N = %d', N));
    end
end

function hileline(src, event)
    % This callback toggles the visibility of the line

    if strcmp(event.Peer.Visible,'on')   % If current line is visible
        event.Peer.Visible = 'off';      %   Set the visibility to 'off'
    else                                 % Else
        event.Peer.Visible = 'on';       %   Set the visibility to 'on'
    end
end