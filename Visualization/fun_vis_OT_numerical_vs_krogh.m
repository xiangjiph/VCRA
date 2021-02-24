function exit_code = fun_vis_OT_numerical_vs_krogh(tmp_cube_str, save_fig_Q)

if nargin < 2
    save_fig_Q = false;
end
persistent DataManager
if isempty(DataManager)
    DataManager = FileManager;
end

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
histogram2(ax_hdl, 'XBinEdges', tmp_cube_str.pO2_dt_hist2.dt_edge, ...
    'YBinEdges', tmp_cube_str.pO2_dt_hist2.pO2_edge, ...
    'BinCounts', tmp_cube_str.pO2_dt_hist2.count, 'DisplayStyle', 'tile');
hold(ax_hdl, 'on');

mean_hdl = plot(ax_hdl, tmp_cube_str.pO2_stat_in_dt_bin.x_bin_val, tmp_cube_str.pO2_stat_in_dt_bin.y_median, 'LineWidth', 4);
ptl005_hdl = plot(ax_hdl, tmp_cube_str.pO2_stat_in_dt_bin.x_bin_val, tmp_cube_str.pO2_stat_in_dt_bin.y_prctile(:, 1), 'LineWidth', 4, ...
    'LineStyle', '-.');
ptl095_hdl = plot(ax_hdl, tmp_cube_str.pO2_stat_in_dt_bin.x_bin_val, tmp_cube_str.pO2_stat_in_dt_bin.y_prctile(:, end), 'LineWidth', 4, ...
    'LineStyle', '-.');
plot_krogh_x = 0 : min(tmp_cube_str.fit_Krogh.d_max, max(tmp_cube_str.local_dt_stat.max));
oxygen_krogh = tmp_cube_str.fit_Krogh.corr_coeff * tmp_cube_str.fit_Krogh.fit_fun_hdl(plot_krogh_x);
plt_krogh = plot(ax_hdl, plot_krogh_x, oxygen_krogh, 'LineWidth', 4, 'LineStyle', '-');

oxygen_krogh_0 = tmp_cube_str.fit_Krogh.fit_fun_hdl(plot_krogh_x);
plt_krogh_0 = plot(ax_hdl, plot_krogh_x, oxygen_krogh_0, 'LineWidth', 4, 'LineStyle', '-');

dot_dt_lm_hdl = scatter(ax_hdl, tmp_cube_str.dt_lm.v, tmp_cube_str.dt_lm.pO2_v, 'y+');
dot_oxy_lm_hdl = scatter(ax_hdl, tmp_cube_str.pO2_lm.dt_v, tmp_cube_str.pO2_lm.v, 'r*');
%
leg_hdl_array = [mean_hdl, ptl005_hdl, ptl095_hdl, plt_krogh, plt_krogh_0, dot_dt_lm_hdl, dot_oxy_lm_hdl];
leg_hdl_string = {'Median', '5% percentile', '95% percentile', ...
    sprintf('Fitting Krogh model:\nCorr. factor: %.2f\nR^2-Adjusted: %.2f', tmp_cube_str.fit_Krogh.corr_coeff, ...
    tmp_cube_str.fit_Krogh.Rsquared.Adjusted), 'Krogh model'...
    'Local d maximum', 'Local u minimum'};
legend(ax_hdl, leg_hdl_array, leg_hdl_string);
ax_hdl.XLabel.String = 'd (\mum)';
ax_hdl.YLabel.String = 'u(d)';
ax_hdl.FontSize = 14;
if save_fig_Q
    vis_folder = fullfile(DataManager.fp_visualization_folder(tmp_cube_str.dataset_name, tmp_cube_str.stack), 'Tissue_pO2_simulation');
    fig_fp = fullfile(vis_folder, sprintf('%s_%s_compare_numerical_to_Krogh_model_240_cube_label_%d.png', ...
        tmp_cube_str.dataset_name, tmp_cube_str.stack, tmp_cube_str.grid_label));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
end
exit_code = 0;
end