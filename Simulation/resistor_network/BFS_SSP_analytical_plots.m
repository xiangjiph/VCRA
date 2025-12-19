clc;clear;close all;
vis_folder = fullfile('.', 'figure', 'BFS_SSP_theory_2025');
fig_name_prefix = 'BFS_SSP';
%% Self-modulation
fig_hdl = figure;
ax_hdl_1 = axes(fig_hdl);
fig_hdl.Position(3:4) = [500, 400];
[eps_ab_array, ratio_array] = meshgrid(-1 : 0.05 : 8, 0 : 0.05 : 1);
flow_ratio_array = (1 + eps_ab_array) ./ (1 + eps_ab_array .* ratio_array) - 1;
[~, contour_hdl] = contourf(ax_hdl_1, eps_ab_array, ratio_array, flow_ratio_array, ...
    [-1 : 0.5 : 0, 0.1 : 0.2 : 0.7, 1, 1.5, 2, 3, 4], 'ShowText', 'on');

ax_hdl_1.XLabel.Interpreter = 'latex';
ax_hdl_1.XLabel.String = '$\hat{\epsilon}_{ij}$';
ax_hdl_1.XTick = [-1, 0 : 2 : 8];

ax_hdl_1.YLabel.Interpreter = 'latex';
ax_hdl_1.YLabel.String = '$f^{\alpha}_{\alpha}$';
ax_hdl_1.FontSize = 12;
grid(ax_hdl_1, 'on');
cbar_hdl = colorbar(ax_hdl_1);
cbar_hdl.Label.Interpreter = 'latex';
cbar_hdl.Label.String = '$\Delta Q_{\alpha}/Q_{\alpha}$';
fig_fp = fullfile(vis_folder, sprintf('%s_self_effect_fabab.png', fig_name_prefix));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Neighbor-modulation flexibility
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
[plt_eps_uAu, plt_uAv2uAu] = meshgrid(-1 : 0.05 : 8, 10.^linspace(-2, log10(1), 100));
plt_z_array = (plt_eps_uAu .* plt_uAv2uAu) ./ (1 + plt_eps_uAu);
plt_z_array(abs(plt_z_array) > 1) = nan;
[~, contour_hdl] = contourf(ax_hdl, plt_eps_uAu, plt_uAv2uAu, plt_z_array);
ax_hdl.YScale = 'log';
cbar_hdl = colorbar(ax_hdl);
cbar_hdl.Label.Interpreter = 'latex';
cbar_hdl.Label.String = '$|\Delta I_{mn} / I_{ij}|$';
cbar_hdl.Label.FontSize = 14;
contour_hdl.ShowText = 'on';
contour_hdl.LevelList = [-3.0 : 1 : -1, -0.1, -0.01, 0, ...
    0.01, 0.05, 0.1 : 0.1 : 0.5];
contour_hdl.LabelSpacing = 300;

ax_hdl.XLabel.Interpreter = 'latex';
ax_hdl.XLabel.String = '$\hat{\epsilon}_{ij}f_{ij}^{ij}$';
ax_hdl.YLabel.Interpreter = 'latex';
ax_hdl.YLabel.String = '$|f^{mn}_{ij} / f^{ij}_{ij}|$';
ax_hdl.FontSize = 12;
grid(ax_hdl, 'on');

fig_fp = fullfile(vis_folder, sprintf('%s_neighbor_effect_fmnab2fabab_epsfabab_cp.png', fig_name_prefix));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%%
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
[eps_ab, fabab] = meshgrid(-1 : 0.05 : 8, 0 : 0.01 : 1);
plt_z_array = eps_ab ./ (1 + eps_ab .* fabab);
[~, contour_hdl] = contourf(ax_hdl, eps_ab, fabab, plt_z_array);

cbar_hdl = colorbar(ax_hdl);
cbar_hdl.Label.Interpreter = 'latex';
cbar_hdl.Label.String = '$\frac{\hat{\epsilon}_{ij}}{1 + \hat{\epsilon}_{ij}f^{ij}_{ij}}$';
cbar_hdl.Label.FontSize = 14;
contour_hdl.ShowText = 'on';
contour_hdl.LevelList = [-4, -2, -1, -0.5, 0, 0.5, 1, 2, 4];
contour_hdl.LabelSpacing = 180;

ax_hdl.XLabel.Interpreter = 'latex';
ax_hdl.XLabel.String = '$\hat{\epsilon}_{ij}$';
ax_hdl.YLabel.Interpreter = 'latex';
ax_hdl.YLabel.String = '$f^{ij}_{ij}$';
ax_hdl.FontSize = 14;
grid(ax_hdl, 'on');

fig_fp = fullfile(vis_folder, sprintf('%s_neighbor_effect_eps_fabab_cp.png', fig_name_prefix));
fun_print_image_in_several_formats(fig_hdl, fig_fp);