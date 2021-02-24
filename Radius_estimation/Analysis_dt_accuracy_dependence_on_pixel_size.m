num_pxl = 30;
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
ax_hdl.XLabel.String = 'Pixel distance';
ax_hdl.YLabel.String = 'Relative error';
ax_hdl.XLim = [0, 15];
ax_hdl.YLim = [1e-3, 1];
ax_hdl.XTick = [0 : 15];
box(ax_hdl, 'on');
grid(ax_hdl, 'on');
ax_hdl.YScale = 'log';
ax_hdl.FontSize = 14;
hold(ax_hdl, 'on');
tmp_mask = strel('disk', num_pxl).Neighborhood;
tmp_mask = padarray(tmp_mask, [1,1], 0);
dist_array = unique(bwdist(~tmp_mask));
plt_2_hdl = plot(ax_hdl , dist_array(2:(end-1)), (dist_array(3:end) - dist_array(1 : (end-2))) ./ (4 * dist_array(2:(end-1))), 'LineWidth', 2);
fig_fp = fullfile(DataManager.fp_visualization_folder('WholeBrain', 'all_stack'), ...
    'paper', 'Pipeline_DT_relative_error_vs_pixel_distance.png');
fun_print_image_in_several_formats(fig_hdl, fig_fp);

%% Estimate using DT
num_pxl = 20;
dist_array = unique(bwdist(~strel('disk', num_pxl).Neighborhood));
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
plot(ax_hdl , dist_array(2:(end-1)), (dist_array(3:end) - dist_array(1 : (end-2))) ./ (4 * dist_array(2:(end-1))), 'LineWidth', 2);
ax_hdl.XLabel.String = 'Pixel distance';
ax_hdl.YLabel.String = 'Relative error';
ax_hdl.XLim = [0, 15];
ax_hdl.YLim = [1e-3, 1];
box(ax_hdl, 'on');
grid(ax_hdl, 'on');
ax_hdl.YScale = 'log';
ax_hdl.FontSize = 14;