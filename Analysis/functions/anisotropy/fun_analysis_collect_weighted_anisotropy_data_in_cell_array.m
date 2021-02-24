function wb_data_str = fun_analysis_collect_weighted_anisotropy_data_in_cell_array(data_cell_array, num_range)

num_valid_cube = numel(data_cell_array);

wb_data_str = struct;
[wb_data_str.num_data, wb_data_str.weight_sum, ...
    wb_data_str.fa, wb_data_str.fa_p, wb_data_str.fa_z, wb_data_str.fa_null_mean, wb_data_str.fa_null_std, ...
    wb_data_str.svd_1, wb_data_str.svd_1_p, wb_data_str.svd_1_z, wb_data_str.svd_1_null_mean, wb_data_str.svd_1_null_std, ...
    ] = deal(nan(num_range, num_valid_cube));
wb_data_str.svd_max_vec = nan(3, num_range, num_valid_cube);
extract_field_prefix = {'fa', 'svd_1'};
extract_field_postfix = {'', '_p', '_z'};
extract_null_postfix = {'mean', 'std'};

%%
num_prefix = numel(extract_field_prefix);
num_postfix = numel(extract_field_postfix);
num_null_field = numel(extract_null_postfix);
start_tic = tic;
for iter_cube = 1 : num_valid_cube
    tmp_data = data_cell_array{iter_cube};
    if ~isempty(tmp_data) && ~isempty(tmp_data.num_data)
        for iter_prefix = 1 : num_prefix
            tmp_prefix_string = extract_field_prefix{iter_prefix};
            for iter_postfix = 1 : num_postfix                 
                tmp_field_name = sprintf('%s%s', tmp_prefix_string, extract_field_postfix{iter_postfix});
                tmp_field_data = tmp_data.(tmp_field_name);
                if ~isempty(tmp_field_data)
                    wb_data_str.(tmp_field_name)(:, iter_cube) = tmp_field_data;
                end
            end
            if isfield(tmp_data, 'null')
               tmp_null_data = tmp_data.null;
               for iter_field = 1 : num_null_field
                   tmp_null_field_data = tmp_null_data.(tmp_prefix_string).(extract_null_postfix{iter_field});
                   if ~isempty(tmp_null_field_data)
                       tmp_new_field_name = sprintf('%s_null_%s', tmp_prefix_string, ...
                           extract_null_postfix{iter_field});
                       wb_data_str.(tmp_new_field_name)(:, iter_cube) = tmp_null_field_data;
                   end
               end
            end
        end
        wb_data_str.svd_max_vec(:, :, iter_cube) = tmp_data.svd_max_vec;
        wb_data_str.num_data(:, iter_cube) = tmp_data.num_data;
        wb_data_str.weight_sum(:, iter_cube) = tmp_data.weight_sum;
    end
    
    if mod(iter_cube, 1000) == 0
        fprintf('Finish processing cube %d / %d (%.2f%%). Elapse time is %.2f seconds.\n', ...
            iter_cube, num_valid_cube, (iter_cube / num_valid_cube * 100), toc(start_tic));
    end
end
%% Transpose the array
field_names = fieldnames(wb_data_str);
for iter_field = 1 : numel(field_names)
    tmp_field_name = field_names{iter_field};
    tmp_field_data = wb_data_str.(tmp_field_name);
    if ismatrix(tmp_field_data)
        wb_data_str.(tmp_field_name) = tmp_field_data.';
    elseif ndims(tmp_field_data) == 3
        wb_data_str.(tmp_field_name) = permute(tmp_field_data, [3, 1, 2]);
    else
        error('Unknown field data type');
    end
end
end