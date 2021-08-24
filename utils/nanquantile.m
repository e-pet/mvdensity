function quant = nanquantile(signal, p, dim)
    assert(ismatrix(signal));
    if (size(signal, 1) == 1) || (nargin > 2 && dim == 2)
        signal = signal';
    end
    n = size(signal, 2);
    quant = zeros(length(p), n);
    for ii = 1:n
        col = signal(:, ii);
        quant(:, ii) = quantile(col(~isnan(col)), p)';
    end
    if nargin > 2 && dim == 2
        quant = quant';
    end
end
