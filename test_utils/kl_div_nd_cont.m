function val = kl_div_nd_cont(pdf1, pdf2, lb, ub)
    % N-Dimensional Kullback-Leibler divergence between pdfs 1 and 2 of continuous random variables.
    % lb and ub should specify the region of integration.
    % For dimensions greater than 3, this method does not work for infinite domains of integration.

    % Notice that for some reason, pdfs in matlab usually expect query points in rows, whereas
    % integration functions expect query points in columns.

    if size(lb, 2) > 1
        lb = lb';
    end
    
    if size(ub, 2) > 1
        ub = ub';
    end
    
    function integrand = kl_integrand(x)
        pdf1x = pdf1(x')';
        integrand = zeros(1, size(x, 2));
        nzmask = pdf1x > 0;
        integrand(nzmask) = pdf1x(nzmask) .* log(pdf1x(nzmask)./pdf2(x(:, nzmask)')');
    end
    
    val = integraln(@kl_integrand, lb, ub);
end