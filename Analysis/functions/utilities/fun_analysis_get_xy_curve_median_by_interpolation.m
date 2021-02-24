function avg_data = fun_analysis_get_xy_curve_median_by_interpolation(x_cell, y_cell)

num_cell = numel(y_cell);
cell_size = size(y_cell);
assert(numel(x_cell) == num_cell);
% Determine the range of x to compute the average value
avg_data.min_x = median(cellfun(@min, x_cell(:)));
avg_data.max_x = median(cellfun(@max, x_cell(:)));

avg_data.interpolate_x_step = median(cellfun(@(x) median(diff(x), 'omitnan'), x_cell(:)), 'omitnan');

avg_data.interpolate_num_step = ceil((avg_data.max_x - avg_data.min_x) / avg_data.interpolate_x_step);
avg_data.interpolate_x = linspace(avg_data.min_x, avg_data.max_x, avg_data.interpolate_num_step);

[avg_data.valid_num_data] = deal(zeros(size(avg_data.interpolate_x)));
tmp_y_data = nan(numel(avg_data.interpolate_x), num_cell);


avg_data.y_interpolation = cell(cell_size);
for iter_cell = 1 : num_cell
    tmp_x = double(x_cell{iter_cell});
    tmp_y = double(y_cell{iter_cell});
    tmp_is_valid_Q = ~isnan(tmp_x) & ~isnan(tmp_y);
    
    avg_data.y_interpolation{iter_cell} = griddedInterpolant(tmp_x(tmp_is_valid_Q), ...
        tmp_y(tmp_is_valid_Q), 'linear', 'none');
    tmp_value = avg_data.y_interpolation{iter_cell}(avg_data.interpolate_x);
    tmp_y_data(:, iter_cell) = tmp_value;
end
avg_data.y_median = median(tmp_y_data, 2, 'omitnan');
end