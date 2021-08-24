addpath('utils');
addpath('test_utils');

logNmax = 6;
iterations_per_N = 10;

rng(1); % set seed

make_dir_if_not_present('./outfiles/');


%% 1-D example
gm = gmdistribution([-1; 1], cat(3, [1], [0.3]), [0.3, 0.7]);
title = '1-D 2-Gaussian mixture example';
export_logN = 3;
analyze_example(gm, title, logNmax, iterations_per_N, './outfiles/MVDensityResults1D.csv', ...
    export_logN, './outfiles/MVDensityExample1D.csv', './outfiles/MVDensityExample1DPDF.csv');


%% 2-D example
gm = gmdistribution([-1, 1; 1, -3], cat(3, [1, 0; 0, 2], [0.3, 0; 0, 0.1]), [0.3, 0.7]);
title = '2-D 2-Gaussian mixture example';
analyze_example(gm, title, logNmax, iterations_per_N, './outfiles/MVDensityResults2D.csv');


%% 5-D example
means = [-1, 1, 3, 0, 0; 0, 0, 2, 3, 4; 2, 2, -2, -3, 3];
covs = cat(3, generateSPDmatrix(5), generateSPDmatrix(5), generateSPDmatrix(5));
weights = [0.3, 0.6, 0.1];
gm = gmdistribution(means, covs, weights);
title = '5-D 3-Gaussian mixture example';
analyze_example(gm, title, logNmax, iterations_per_N, './outfiles/MVDensityResults5D.csv');