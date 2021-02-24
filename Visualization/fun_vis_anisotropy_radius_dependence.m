function fig_hdl = fun_vis_anisotropy_radius_dependence(anisotropy_all_vw_range, vis_range)
%% Visualization
if nargin < 2
    vis_range = [0, inf];
end

fa_p = max(5e-5, anisotropy_all_vw_range.fa_p);
min2max_p = max(5e-5, anisotropy_all_vw_range.svd_min2max_p);
svd_1_p = max(5e-5, anisotropy_all_vw_range.svd_1_p);

radius_selection_max = anisotropy_all_vw_range.select_r_max;
radius_selection_min = anisotropy_all_vw_range.select_r_min;
plot_list_ind = find(radius_selection_min == vis_range(1) & ...
    radius_selection_max <= vis_range(2));
% plot_list_ind = 1 : numel(radius_selection_max);
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [1, 1.6];
ax_hdl_1 = subplot(4, 1, 1);
plot(ax_hdl_1, radius_selection_max(plot_list_ind), ...
    anisotropy_all_vw_range.weight_sum(plot_list_ind) ./...
    max(anisotropy_all_vw_range.weight_sum(plot_list_ind)));
ax_hdl_1.YLim = [0, 1];
ax_hdl_1.YLabel.String = sprintf('Fraction of total %s', ...
    anisotropy_all_vw_range.weight_method);
ax_hdl_1.XAxis.Visible = 'off';
ax_hdl_1.XScale = 'log';
grid(ax_hdl_1, 'on');

ax_hdl_2 = subplot(4, 1, 2);
yyaxis(ax_hdl_2, 'left');
tmp_plt_m = plot(ax_hdl_2, radius_selection_max(plot_list_ind), ...
    anisotropy_all_vw_range.fa(plot_list_ind));
hold(ax_hdl_2, 'on');
tmp_plt_n = errorbar(ax_hdl_2, radius_selection_max(plot_list_ind), ...
    anisotropy_all_vw_range.null.fa.mean(plot_list_ind), ...
    anisotropy_all_vw_range.null.fa.std(plot_list_ind));
ax_hdl_2.YLim = [0, 1];
ax_hdl_2.YLabel.String = 'FA';
yyaxis(ax_hdl_2, 'right');
plot(ax_hdl_2, radius_selection_max(plot_list_ind), ...
    fa_p(plot_list_ind));
ax_hdl_2.YLabel.String = 'FA_p';
ax_hdl_2.YScale = 'log';
legend(ax_hdl_2, [tmp_plt_m, tmp_plt_n], 'Measured', 'Null', 'Location', 'northeast');
ax_hdl_2.XAxis.Visible = 'off';
ax_hdl_2.XScale = 'log';
grid(ax_hdl_2, 'on');

ax_hdl_3 = subplot(4, 1, 3);
yyaxis(ax_hdl_3, 'left');
tmp_plt_m = plot(ax_hdl_3, radius_selection_max(plot_list_ind), ...
    anisotropy_all_vw_range.svd_1(plot_list_ind));
hold(ax_hdl_3, 'on');
tmp_plt_n = errorbar(ax_hdl_3, radius_selection_max(plot_list_ind), ...
    anisotropy_all_vw_range.null.svd_1.mean(plot_list_ind), ...
    anisotropy_all_vw_range.null.svd_1.std(plot_list_ind));
ax_hdl_3.YLim = [0, 1];
ax_hdl_3.YLabel.String = 'PCV1';
yyaxis(ax_hdl_3, 'right');
plot(ax_hdl_3, radius_selection_max(plot_list_ind), ...
    svd_1_p(plot_list_ind));
ax_hdl_3.YLabel.String = 'PCV1_p';
ax_hdl_3.YScale = 'log';
ax_hdl_3.XAxis.Visible = 'off';
legend(ax_hdl_3, [tmp_plt_m, tmp_plt_n], 'Measured', 'Null', 'Location', 'northeast');
ax_hdl_3.XScale = 'log';
grid(ax_hdl_3, 'on');

ax_hdl_4 = subplot(4, 1, 4);
yyaxis(ax_hdl_4, 'left');
tmp_plt_m = plot(ax_hdl_4, radius_selection_max(plot_list_ind), ...
    anisotropy_all_vw_range.svd_min2max(plot_list_ind));
hold(ax_hdl_4, 'on');
tmp_plt_n = errorbar(ax_hdl_4, radius_selection_max(plot_list_ind), ...
    anisotropy_all_vw_range.null.svd_min2max.mean(plot_list_ind), ...
    anisotropy_all_vw_range.null.svd_min2max.std(plot_list_ind));
ax_hdl_4.YLim = [0, 1];
ax_hdl_4.YLabel.String = 'PCV min2max';
yyaxis(ax_hdl_4, 'right');
plot(ax_hdl_4, radius_selection_max(plot_list_ind), ...
    min2max_p(plot_list_ind));
ax_hdl_4.YScale = 'log';
ax_hdl_4.YLabel.String = 'PCV min2max_p';
ax_hdl_4.XLabel.String = 'Capillary radius cutoff (\mum)';
ax_hdl_4.Box = 'off';
legend(ax_hdl_4, [tmp_plt_m, tmp_plt_n], 'Measured', 'Null', 'Location', 'northeast');
ax_hdl_4.XScale = 'log';
grid(ax_hdl_4, 'on');
set(fig_hdl.Children, 'FontSize', 12);
set(fig_hdl.Children, 'LineWidth', 1);

%     fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), ...
%         'Anisotropy_r_dependence', sprintf('%s_%s_%s_%d_%d_%d_%s_weighted_anisotropy_r_dependence.png', ...
%         dataset_name, stack, skel_version, tmp_grid_sub, anisotropy_test_method));
%     fun_print_image_in_several_formats(fig_hdl, fig_fp);

%% Visualize the consistence of the principle eigenvector
cos_pcvl = real(acosd(abs(fun_pcos(anisotropy_all_vw_range.svd_max_vec'))));
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
imagesc(ax_hdl, cos_pcvl);
cbar_hdl = colorbar(ax_hdl);
ax_hdl.CLim = [0, 90];
ax_hdl.ColorScale = 'linear';
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.XTickLabel = arrayfun(@(x) num2str(x, '%.1f'), radius_selection_max(ax_hdl.XTick), 'UniformOutput', false);
ax_hdl.YTickLabel = arrayfun(@(x) num2str(x, '%.1f'), radius_selection_max(ax_hdl.YTick), 'UniformOutput', false);
end