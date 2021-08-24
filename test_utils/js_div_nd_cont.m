function val = js_div_nd_cont(pdf1, pdf2, lb, ub)
    % N-Dimensional Jensen-Shannon divergence between pdfs 1 and 2 of continuous random variables.
    % lb and ub should specify the region of integration.
    % For dimensions greater than 3, this method does not work for infinite domains of integration.
    m = @(x) (pdf1(x)+pdf2(x)) ./ 2;
    val = 0.5 * kl_div_nd_cont(pdf1, m, lb, ub) + 0.5 * kl_div_nd_cont(pdf2, m, lb, ub);
end