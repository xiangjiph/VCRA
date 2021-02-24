function avg_data = fun_analysis_get_xy_curve_avgNstd_by_interpolation(x_cell, y_cell, int_x_list, extrapolation_method)

num_cell = numel(y_cell);
cell_size = size(y_cell);
assert(numel(x_cell) == num_cell);
%% Determine the range of x to compute the average value
if nargin < 3    
    avg_data.min_x = median(cellfun(@min, x_cell(:)));
    avg_data.max_x = median(cellfun(@max, x_cell(:)));
    
    avg_data.interpolate_x_step = abs(median(cellfun(@(x) median(diff(x), 'omitnan'), x_cell(:)), 'omitnan'));
    
    avg_data.interpolate_num_step = ceil((avg_data.max_x - avg_data.min_x) / avg_data.interpolate_x_step) + 1;
    avg_data.interpolate_x = linspace(avg_data.min_x, avg_data.max_x, avg_data.interpolate_num_step);
else
    avg_data.interpolate_x = int_x_list;
end
if nargin < 4
    extrapolation_method = 'none';
end
%%
[avg_data.valid_num_data, avg_data.y_squared_avg, ...
    avg_data.y_avg] = deal(zeros(size(avg_data.interpolate_x)));

avg_data.y_interpolation = cell(cell_size);
for iter_cell = 1 : num_cell
    tmp_x = double(x_cell{iter_cell});
    tmp_y = double(y_cell{iter_cell});
    is_valid_Q = ~isnan(tmp_x) & ~isnan(tmp_y);
    tmp_x = tmp_x(is_valid_Q);
    tmp_y = tmp_y(is_valid_Q);    
    if ~issorted(tmp_x, 'ascend')
        [tmp_x, tmp_sort_idx] = sort(tmp_x, 'ascend');
        tmp_y = tmp_y(tmp_sort_idx);
    end
    [tmp_x, tmp_unique_idx] = unique(tmp_x, 'stable');
    tmp_y = tmp_y(tmp_unique_idx);
    
    avg_data.y_interpolation{iter_cell} = griddedInterpolant(tmp_x, tmp_y, 'linear', extrapolation_method);
    tmp_value = avg_data.y_interpolation{iter_cell}(avg_data.interpolate_x);
    tmp_is_finite_Q = isfinite(tmp_value);
    avg_data.valid_num_data = avg_data.valid_num_data + tmp_is_finite_Q;
    tmp_value(~tmp_is_finite_Q) = 0;
    avg_data.y_avg = avg_data.y_avg + tmp_value;
    avg_data.y_squared_avg = avg_data.y_squared_avg + tmp_value .^ 2;
end
avg_data.y_avg = avg_data.y_avg ./ avg_data.valid_num_data;
avg_data.y_squared_avg = avg_data.y_squared_avg ./ avg_data.valid_num_data;

avg_data.y_std = sqrt(avg_data.y_squared_avg - avg_data.y_avg .^ 2);

avg_data.y_avg_interpolation = griddedInterpolant(avg_data.interpolate_x, avg_data.y_avg, 'linear', extrapolation_method);
avg_data.y_std_interpolation = griddedInterpolant(avg_data.interpolate_x, avg_data.y_std, 'linear', extrapolation_method);
end