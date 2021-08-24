function [bin_flags, bin_counts, bin_flag_idx, data_flag] = hist_bin_counts(data, bin_lims, varargin)

    [N, m] = size(data);

    assert(m >= 1)
    assert(length(varargin) == m - 1)
    
    bin_lims_cl = cell(1, m);
    bin_lims_cl{1} = bin_lims;
    
    for ii = 2:m
        bin_lims_cl{ii} = varargin{ii-1};
    end
    cellfun(@(arr) assert(issorted(arr)), bin_lims_cl);
    
    nan_rows = any(isnan(data), 2);
    
    % data_flag has the same dimension as data and contains the
    % corresponding flags for each dimension of each sample.
    data_flag = zeros(size(data));
    for ii = 1:m
        for jj=1:length(bin_lims_cl{ii})
            data_flag(:, ii) = data_flag(:, ii) + (data(:, ii) > bin_lims_cl{ii}(jj));
        end
        assert(all(data_flag(~nan_rows, ii) >= 0));
        assert(all(data_flag(~nan_rows, ii) <= length(bin_lims_cl{ii})));
    end
    
    assert(all(size(data_flag) == size(data)));
    data_flag(nan_rows, :) = nan;
    
    % Now count the occurences of each combination of flags, i.e. (for three
    % dimensions), the number of data samples with flag combination "0-2-1" 
    % and so on.
    % bin_flags is a matrix that contains unique flag
    % combination i in row i.
    % bin_flag_idx contains in row i the index of the unique flag
    % combination in bin_flags that is observed in row i of
    % data_flag.
    bin_flag_idx = nan * ones(size(data, 1), 1);
    [bin_flags, ~, bin_flag_idx(~nan_rows)] = unique(data_flag(~nan_rows, :), 'rows');
    [~, nf] = size(bin_flags);
    assert(nf == m);
    % Vector that contains the number of occurences of flag combination i
    % in position i
    bin_counts = accumarray(bin_flag_idx(~nan_rows), 1);
    assert(sum(bin_counts) == N - sum(nan_rows))
end