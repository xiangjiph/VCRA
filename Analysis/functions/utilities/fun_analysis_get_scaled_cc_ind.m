function cc_ind = fun_analysis_get_scaled_cc_ind(cc_ind, mask_size, scale, target_size)

if isempty(cc_ind)
    return;
else
    num_cc = numel(cc_ind);
    if nargin < 4
        target_size = round(mask_size .* scale);
    end
    assert(ndims(target_size) == ndims(mask_size));
    for iter_cc = 1 : num_cc
        tmp_sub = fun_ind2sub(mask_size, cc_ind{iter_cc});
        tmp_sub = round(tmp_sub .* scale);
        tmp_sub = min(max(1, tmp_sub), target_size);
        cc_ind{iter_cc} = sub2ind(target_size, tmp_sub(:, 1), tmp_sub(:, 2), ...
            tmp_sub(:, 3));
    end
end
end