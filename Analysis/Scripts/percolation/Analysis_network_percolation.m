set_env;
clc;clear;close all;
vessel_graph = load('test_graph.mat');
vessel_graph_wo_ep = fun_graph_pruning_internal_short_hairs(vessel_graph, inf, 0);
graph_ud_str = fun_analysis_get_connectivity_graph(vessel_graph_wo_ep);
graph_ud = graph_ud_str.graph_w;

graph_degree = degree(graph_ud);
phi_c = (mean(graph_degree)) / (mean(graph_degree.^2) - mean(graph_degree));
%%
percolation_p =  0.05 : 0.1 : 0.95;
percolation_str = fun_analysis_percolation_transition_by_bond_removal(graph_ud, percolation_p,  ...
    50, false);
%% Percolation transition 
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3,1];
ax_hdl = subplot(1,3,1);
plot(ax_hdl, percolation_str.bond_occupancy_p , percolation_str.avg.largest_cc_fraction);
ax_hdl.YLabel.String = 'Largest cluster fraction (s)';
ax_hdl.XLabel.String = 'Occupation probability (p)';
ax_hdl.XLim = [0,1];
% Find the value of remove fraction s.t. largest cc ratio = 0.5 by linear
% interpolation 
[int_x, tmp_ind] = sort(percolation_str.avg.largest_cc_fraction, 'ascend');
int_y = percolation_str.bond_occupancy_p(tmp_ind);
cc_ratio_2_p = griddedInterpolant(int_x, int_y);

ratio_half_th = cc_ratio_2_p(0.5);
ratio_095_th = cc_ratio_2_p(0.95);
ratio_005_th = cc_ratio_2_p(0.05);
ratio_001_th = cc_ratio_2_p(0.01);
legend(ax_hdl, sprintf('p(s = 0.05) = %.4f\np(s = 0.50) = %.4f\np(s = 0.95) = %.4f', ...
    ratio_005_th, ratio_half_th, ratio_095_th), 'Location', 'northwest');
%% Divergent exponent near threshold 
% ax_hdl_2 = axes(figure);
ax_hdl_2 = subplot(1,3,2);
s_range = [0.025, 0.5];
p_range = cc_ratio_2_p(s_range);

loglog_fit_Q = (percolation_p <= p_range(2) & percolation_p >= p_range(1));
fit_x = percolation_p(loglog_fit_Q) - p_range(1);
fit_y = log10(avg_largest_cc_ratio(loglog_fit_Q));
scatter(ax_hdl_2, fit_x, fit_y);
lin_log_fit = fitlm(fit_x, fit_y);
hold(ax_hdl_2, 'on');
plt_hdl = plot(ax_hdl_2, fit_x, fit_x * lin_log_fit.Coefficients.Estimate(2) + lin_log_fit.Coefficients.Estimate(1));
% ax_hdl_2.YLim = [0,1];
ax_hdl_2.YLabel.String = 'log_{10}(s)';
ax_hdl_2.XLabel.String = sprintf('p - p_c(s = %.3f)', s_range(1));
legend(plt_hdl, sprintf('p_c: %.4f\nSlope: %.2f \\pm %.2f\nIntercept: %.2f \\pm %.2f\nR-squared: %.3f\nData size: %d', ...
    p_range(1), ...
    lin_log_fit.Coefficients.Estimate(2), lin_log_fit.Coefficients.SE(2), ...
    lin_log_fit.Coefficients.Estimate(1), lin_log_fit.Coefficients.SE(1), ...
    lin_log_fit.Rsquared.Adjusted, lin_log_fit.NumObservations), 'Location', 'northwest');
%% 
ax_hdl_3 = subplot(1,3,3);
scatter(ax_hdl_3, percolation_p(loglog_fit_Q) - p_range(1), ...
    avg_largest_cc_ratio(loglog_fit_Q));
ax_hdl_3.XScale = 'log';
ax_hdl_3.YScale = 'log';
ax_hdl_3.XLabel.String = sprintf('p - p_c(s = %.3f)', s_range(1));
ax_hdl_3.YLabel.String = 's';
% grid(ax_hdl_3, 'on');
%%
DataManager = FileManager;
dataset_name = vessel_graph.info.dataset_name;
stack = vessel_graph.info.stack;
grid_version = vessel_graph.info.grid_version;
% fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'percolation', ...
%     sprintf('%s_%s_%s_%d_giant_cc_vs_link_rm_fraction.png', dataset_name, stack, grid_version, vessel_graph.info.combined_grid_xyz_label));
fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'percolation', ...
    sprintf('%s_%s_%s_%d_giant_cc_node_size_vs_p_log_lin_fit.png', dataset_name, stack, grid_version, vessel_graph.info.combined_grid_xyz_label));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
