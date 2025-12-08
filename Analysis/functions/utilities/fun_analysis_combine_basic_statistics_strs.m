function combined_stat_str = fun_analysis_combine_basic_statistics_strs(stat_str_array)
%

scalar_field_name = {'num_data', 'min', 'max', 'mean', 'median', 'sum', 'std', ...
    'range', 'hist_num_edge'};
array_field_name = setdiff(fieldnames(stat_str_array), scalar_field_name);
% array_field_name = {'hist_count', 'hist_edge', 'hist_bin_val', ...
%     'hist_probability', 'hist_pdf', 'hist_bin_size', ...
%     'hist_cumulate_count', 'hist_cdf', 'prctile_th', 'prctile_val'};
num_strs = numel(stat_str_array);
combined_stat_str = struct;
for iter_field = 1 : numel(scalar_field_name)
    tmp_field_name = scalar_field_name{iter_field};
    tmp_vec = nan(1, num_strs);
    for iter_col = 1 : num_strs
        tmp_scalar = stat_str_array(iter_col).(tmp_field_name);
        if ~isempty(tmp_scalar)
           tmp_vec(iter_col) = tmp_scalar;
        end
    end
    combined_stat_str.(tmp_field_name) = tmp_vec;
end

for iter_field = 1 : numel(array_field_name)
    tmp_field_name = array_field_name{iter_field};
    combined_stat_str.(tmp_field_name) = {stat_str_array.(tmp_field_name)};
end
combined_stat_str.prctile_th = combined_stat_str.prctile_th{1};
combined_stat_str.prctile_val = cat(1, combined_stat_str.prctile_val{:});

end