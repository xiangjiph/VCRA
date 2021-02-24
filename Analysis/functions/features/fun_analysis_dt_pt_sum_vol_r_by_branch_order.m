function order_vol_r = fun_analysis_dt_pt_sum_vol_r_by_branch_order(branch_order, volume_ratio, order_bin_edge)

num_data = numel(volume_ratio);
num_bin = numel(order_bin_edge) - 1;

if num_data > 0 && num_data == numel(branch_order)
    list_idx_cell = fun_bin_data_to_idx_list_by_edges(branch_order, order_bin_edge, true);
    order_vol_r = nan(num_bin, 1);
    assert(numel(list_idx_cell) == num_bin, 'Mismatch array size');
    for iter_bin = 1 : num_bin
       order_vol_r(iter_bin) = sum(volume_ratio(list_idx_cell{iter_bin})); 
    end    
else
    order_vol_r = nan(num_bin, 1);
end
end