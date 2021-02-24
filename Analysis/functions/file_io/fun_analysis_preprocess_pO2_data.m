function pO2_stat = fun_analysis_preprocess_pO2_data(wholebrain_pO2_data, vector_num_col)

if nargin < 2
    vector_num_col = 1;
end

pO2_stat = struct;
num_cubes = numel(wholebrain_pO2_data.dt_lm);
is_empty_cell_Q = wholebrain_pO2_data.has_no_data_Q;


% pO2_stat.fit_krogh.coefficient = wholebrain_pO2_data.fit_krogh_coeff;

extract_field_name = cell(3, 0);
extract_field_name(:, end+1) = {'dt_lm', 'dt_mean', 'vector'};
extract_field_name(:, end+1) = {'dt_lm', 'dt_median', 'vector'};
extract_field_name(:, end+1) = {'dt_lm', 'pO2_mean', 'vector'};
extract_field_name(:, end+1) = {'dt_lm', 'pO2_median', 'vector'};
extract_field_name(:, end+1) = {'dt_lm', 'dt_mean_prctile', 'vector'};
extract_field_name(:, end+1) = {'dt_lm', 'dt_median_prctile', 'vector'};
extract_field_name(:, end+1) = {'dt_lm', 'pO2_mean_prctile', 'vector'};
extract_field_name(:, end+1) = {'dt_lm', 'pO2_median_prctile', 'vector'};

extract_field_name(:, end+1) = {'pO2_lm', 'dt_mean', 'vector'};
extract_field_name(:, end+1) = {'pO2_lm', 'dt_median', 'vector'};
extract_field_name(:, end+1) = {'pO2_lm', 'pO2_mean', 'vector'};
extract_field_name(:, end+1) = {'pO2_lm', 'pO2_median', 'vector'};
extract_field_name(:, end+1) = {'pO2_lm', 'dt_mean_prctile', 'vector'};
extract_field_name(:, end+1) = {'pO2_lm', 'dt_median_prctile', 'vector'};
extract_field_name(:, end+1) = {'pO2_lm', 'pO2_mean_prctile', 'vector'};
extract_field_name(:, end+1) = {'pO2_lm', 'pO2_median_prctile', 'vector'};

% extract_field_name(:, end+1) = {'paired_extrema', 'dist', 'cell'};
% extract_field_name(:, end+1) = {'paired_extrema', 'pO2', 'cell'};
% extract_field_name(:, end+1) = {'paired_extrema', 'dt', 'cell'};

extract_field_name(:, end+1) = {'pO2_stat_in_dt_bin', 'y_median', 'cell'};
extract_field_name(:, end+1) = {'pO2_stat_in_dt_bin', 'y_mean', 'cell'};
extract_field_name(:, end+1) = {'pO2_stat_in_dt_bin', 'y_std', 'cell'};
extract_field_name(:, end+1) = {'pO2_stat_in_dt_bin', 'x_bin_val', 'cell'};

extract_field_name(:, end+1) = {'local_dt_stat', 'mean', 'scalar'};
extract_field_name(:, end+1) = {'local_dt_stat', 'median', 'scalar'};
extract_field_name(:, end+1) = {'local_dt_stat', 'std', 'scalar'};
extract_field_name(:, end+1) = {'local_dt_stat', 'hist_cdf', 'cell'};
extract_field_name(:, end+1) = {'local_dt_stat', 'hist_pdf', 'cell'};
extract_field_name(:, end+1) = {'local_dt_stat', 'hist_edge', 'cell'};
extract_field_name(:, end+1) = {'local_dt_stat', 'hist_bin_val', 'cell'};
extract_field_name(:, end+1) = {'local_dt_stat', 'prtl2val_itp', 'cell'};
extract_field_name(:, end+1) = {'local_dt_stat', 'val2ptrl_itp', 'cell'};


extract_field_name(:, end+1) = {'local_pO2_stat', 'mean', 'scalar'};
extract_field_name(:, end+1) = {'local_pO2_stat', 'median', 'scalar'};
extract_field_name(:, end+1) = {'local_pO2_stat', 'std', 'scalar'};
extract_field_name(:, end+1) = {'local_pO2_stat', 'hist_cdf', 'cell'};
extract_field_name(:, end+1) = {'local_pO2_stat', 'hist_pdf', 'cell'};
extract_field_name(:, end+1) = {'local_pO2_stat', 'hist_edge', 'cell'};
extract_field_name(:, end+1) = {'local_pO2_stat', 'hist_bin_val', 'cell'};
extract_field_name(:, end+1) = {'local_pO2_stat', 'prtl2val_itp', 'cell'};
extract_field_name(:, end+1) = {'local_pO2_stat', 'val2ptrl_itp', 'cell'};
num_extract_field = size(extract_field_name, 2);
for iter_field = 1 : num_extract_field
    tmp_fn_l1 = extract_field_name{1, iter_field};
    tmp_fn_l2 = extract_field_name{2, iter_field};
    tmp_data_type = extract_field_name{3, iter_field};
    
    switch tmp_data_type
        case 'cell'
            pO2_stat.(tmp_fn_l1).(tmp_fn_l2) = cell(num_cubes, 1);
            pO2_stat.(tmp_fn_l1).(tmp_fn_l2)(~is_empty_cell_Q) = cellfun(@(x) x.(tmp_fn_l2), wholebrain_pO2_data.(tmp_fn_l1)(~is_empty_cell_Q), 'UniformOutput', false);
        case 'scalar'
            pO2_stat.(tmp_fn_l1).(tmp_fn_l2) = nan(num_cubes, 1);
            pO2_stat.(tmp_fn_l1).(tmp_fn_l2)(~is_empty_cell_Q) = cellfun(@(x) x.(tmp_fn_l2), wholebrain_pO2_data.(tmp_fn_l1)(~is_empty_cell_Q), 'UniformOutput', true);
        case 'vector'
            tmp_mat = nan(num_cubes, vector_num_col);
            tmp_cell = cellfun(@(x) x.(tmp_fn_l2)', wholebrain_pO2_data.(tmp_fn_l1)(~is_empty_cell_Q), 'UniformOutput', false);
            tmp_mat(~is_empty_cell_Q, :) = cat(1, tmp_cell{:});
            pO2_stat.(tmp_fn_l1).(tmp_fn_l2) = tmp_mat;
    end
end

end