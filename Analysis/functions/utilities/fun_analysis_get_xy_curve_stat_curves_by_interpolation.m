function avg_data = fun_analysis_get_xy_curve_stat_curves_by_interpolation(x_cell, y_cell, int_x_list, etp_method)

num_cell = numel(y_cell);
cell_size = size(y_cell);
assert(numel(x_cell) == num_cell, 'Number of cells in x is different from y');
%% Determine the range of x to compute the average value
if nargin < 3 || isempty(int_x_list)
    avg_data.min_x = median(cellfun(@min, x_cell(:)));
    avg_data.max_x = median(cellfun(@max, x_cell(:)));
    
    avg_data.interpolate_x_step = median(cellfun(@(x) median(diff(x), 'omitnan'), x_cell(:)), 'omitnan');
    
    avg_data.interpolate_num_step = ceil((avg_data.max_x - avg_data.min_x) / avg_data.interpolate_x_step) + 1;
    avg_data.interpolate_x = linspace(avg_data.min_x, avg_data.max_x, avg_data.interpolate_num_step);
else
    if iscolumn(int_x_list)
        int_x_list = int_x_list.';
    end
    avg_data.interpolate_x = int_x_list;
    avg_data.interpolate_num_step = numel(avg_data.interpolate_x);
end
if nargin < 4
    etp_method = 'none';
end
%%
[avg_data.valid_num_data, avg_data.y_squared_avg, ...
    avg_data.y_avg] = deal(zeros(size(avg_data.interpolate_x)));

avg_data.y_interpolation = cell(cell_size);
avg_data.y_interpolation_value = cell(cell_size);
for iter_cell = 1 : num_cell
    tmp_x = double(x_cell{iter_cell});
    tmp_y = double(y_cell{iter_cell});
    is_valid_Q = ~isnan(tmp_x) & ~isnan(tmp_y);
    tmp_x = tmp_x(is_valid_Q);
    tmp_y = tmp_y(is_valid_Q);    
    [tmp_x, tmp_ind] = unique(tmp_x, 'sorted');
    tmp_y = tmp_y(tmp_ind);
    avg_data.y_interpolation{iter_cell} = griddedInterpolant(tmp_x, tmp_y, 'linear', etp_method);
    tmp_value = avg_data.y_interpolation{iter_cell}(avg_data.interpolate_x);
    avg_data.y_interpolation_value{iter_cell} = tmp_value;
    tmp_is_finite_Q = isfinite(tmp_value);
    avg_data.valid_num_data = avg_data.valid_num_data + tmp_is_finite_Q;
    tmp_value(~tmp_is_finite_Q) = 0;
    avg_data.y_avg = avg_data.y_avg + tmp_value;
    avg_data.y_squared_avg = avg_data.y_squared_avg + tmp_value .^ 2;
end
avg_data.y_avg = avg_data.y_avg ./ avg_data.valid_num_data;
avg_data.y_squared_avg = avg_data.y_squared_avg ./ avg_data.valid_num_data;

avg_data.y_std = sqrt(avg_data.y_squared_avg - avg_data.y_avg .^ 2);

avg_data.y_avg_interpolation = griddedInterpolant(avg_data.interpolate_x, avg_data.y_avg, 'linear', etp_method);
avg_data.y_std_interpolation = griddedInterpolant(avg_data.interpolate_x, avg_data.y_std, 'linear', etp_method);
%%
itp_y_val_mat = cat(1, avg_data.y_interpolation_value{:});
prctile_val = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95];
num_prctile = numel(prctile_val);
prctile_value_array = nan(num_prctile, avg_data.interpolate_num_step);
for iter_step = 1 : avg_data.interpolate_num_step
   tmp_val = itp_y_val_mat(:, iter_step);
   tmp_val = sort(tmp_val(~isnan(tmp_val)), 'ascend');
   tmp_num_elem = numel(tmp_val);
   if tmp_num_elem > 0
       tmp_prctile_ind = ceil(tmp_num_elem .* prctile_val);
       prctile_value_array(:, iter_step) = tmp_val(tmp_prctile_ind);
   end
end
avg_data.prctile = prctile_val;
avg_data.prctile_value = prctile_value_array';
avg_data.prctile_interpolation = cell(num_prctile, 1);
for iter_ptl = 1 : num_prctile
    tmp_int_y = avg_data.prctile_value(:, iter_ptl);
    tmp_int_x = avg_data.interpolate_x;
    tmp_valid_Q = isfinite(tmp_int_y(:)) & isfinite(tmp_int_x(:));
    tmp_int_x = tmp_int_x(tmp_valid_Q);
    tmp_int_y = tmp_int_y(tmp_valid_Q);
    avg_data.prctile_interpolation{iter_ptl} = griddedInterpolant(...
        tmp_int_x, tmp_int_y, 'linear', etp_method);
end
end