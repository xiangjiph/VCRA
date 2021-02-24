function [mean_vec, median_vec] = fun_analysis_get_mean_N_median_in_each_cell(data_cell)

num_cell = numel(data_cell);
[mean_vec, median_vec] = deal(nan(num_cell, 1));
for iter_cell = 1 : num_cell
    tmp_data = data_cell{iter_cell};
    is_valid_Q = isfinite(tmp_data);
    mean_vec(iter_cell) = mean(tmp_data(is_valid_Q));
    median_vec(iter_cell) = median(tmp_data(is_valid_Q));    
end


end