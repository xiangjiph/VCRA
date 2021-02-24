function cc_median_r = fun_analysis_get_cc_median_radius(cc_ind, ind_2_r_map)

assert(iscell(cc_ind));
num_cc = numel(cc_ind);
cc_median_r = nan(num_cc, 1);
is_sparse_Q = issparse(ind_2_r_map);
for iter_cc = 1 : num_cc
    tmp_r = ind_2_r_map(cc_ind{iter_cc});
    if is_sparse_Q
        tmp_r = full(tmp_r);
    end
    cc_median_r(iter_cc) = median(tmp_r(tmp_r > 0));
end
end