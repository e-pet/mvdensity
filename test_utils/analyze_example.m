function analyze_example(test_pdf, title, logNmax, iterations_per_N, export_file, export_logN, ...
                         export_example_file_samples, export_example_file_pdf)

if nargin < 6
    export_logN = NaN;
end

pdf_est_method = @(xn, do_plot) est_pdf_grid(xn, [], [], [], do_plot);

disp('');
fprintf('\n----- %s -----\n', title);

logN = 1:logNmax;
jsds_log_med = zeros(length(logN), 1);
jsds_log_min = zeros(length(logN), 1);
jsds_log_max = zeros(length(logN), 1);
for ii = 1:length(logN)
    N = 10^logN(ii);
    jsds = zeros(iterations_per_N, 1);
    for jj = 1:iterations_per_N
        if logN(ii) == export_logN && jj == 1
            jsds(jj) = run_example(test_pdf, pdf_est_method, N, jj==1, title, ...
                export_example_file_samples, export_example_file_pdf);
        else
            jsds(jj) = run_example(test_pdf, pdf_est_method, N, jj==1, title);
        end
    end
    jsds_log_med(ii) = median(log(jsds));
    jsds_log_min(ii) = min(log(jsds));
    jsds_log_max(ii) = max(log(jsds));
end

fig = figure('NumberTitle', 'off', 'Name', sprintf('JSD: %s', title));
errorbar(logN, jsds_log_med, jsds_log_med-jsds_log_min, jsds_log_max-jsds_log_med, 'LineWidth', 1);

xlabel('log10(N)');
ylabel('log10(JSD)');
legend('MAKIMA Estimate');
gridTitle(fig, sprintf('%s, pdf JSD', title));

if nargin > 4
    writematrix([10.^(logN') exp(jsds_log_med) exp(jsds_log_med)-exp(jsds_log_min) exp(jsds_log_max)-exp(jsds_log_med)], ...
        export_file);
end