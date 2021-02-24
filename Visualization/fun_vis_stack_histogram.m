function [fig_hdl, ax_hdl, varargout] = fun_vis_stack_histogram(data_x, data_y, bin_x_edge, bin_y_edge)

tmp_not_nan_Q = ~isnan(data_x) & ~isnan(data_y);
data_x = data_x(tmp_not_nan_Q);
data_y = data_y(tmp_not_nan_Q);

[bin_x_ind, ~] = fun_bin_data_to_idx_list_by_edges(data_x, bin_x_edge, true);

tmp_bin_cell = cellfun(@(x) histcounts(data_y(x), bin_y_edge, 'Normalization', 'count'), ...
    bin_x_ind, 'UniformOutput', false);
tmp_bin_cell = cat(1, tmp_bin_cell{:});
% Normalization
tmp_bin_cell = tmp_bin_cell ./ numel(data_x);
assert((sum(tmp_bin_cell(:)) - 1) <= eps);

num_edge = numel(bin_x_edge);
tmp_bin_val = (1 : num_edge - 1) - 0.5;

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
bar_hdl = bar(ax_hdl, tmp_bin_val, tmp_bin_cell, 'stacked');
[bar_hdl.BarWidth] = deal(1);
ax_hdl.XTick = 0 : num_edge;
ax_hdl.XTickLabel = arrayfun(@(x) num2str(x, '%.2f'), bin_x_edge, 'UniformOutput', false);

if nargout > 2
    stat_str = struct;
    stat_str.dim_1_bin_edge = bin_x_edge;
    stat_str.dim_2_bin_edge = bin_y_edge;
    stat_str.probability_mat = tmp_bin_cell;
    varargout{1} = stat_str;
end
end