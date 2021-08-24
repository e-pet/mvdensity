A simple method for estimating multivariate probability densities based on sample data. The method works efficiently in multiple dimensions, with many datapoints (tested up to a million), and is also cheap to *evaluate* on many datapoints. It works simply by computing a multivariate histogram and then smoothing or interpolating the histogram counts.

For more details, see
- the main function, est_pdf_grid.m
- the doc file documentation.pdf
- the test examples implemented in run_pdf_examples.m

-- Eike Petersen, August 2021