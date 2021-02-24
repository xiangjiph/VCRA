function [fig_hdl, ax_hdl]= fun_analysis_radius_calibration_vis(calibration_str, stat_name)

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
histogram2(ax_hdl, calibration_str.Moving_r_um, calibration_str.Fixed_r_um, ...
    calibration_str.Hist_edge, calibration_str.Hist_edge, 'DisplayStyle', 'tile');
hold(ax_hdl, 'on');
switch stat_name
    case {'Median', 'median'}
        errorbar(ax_hdl, calibration_str.Fixed_binned_by_Moving.x_bin_val, ...
            calibration_str.Fixed_binned_by_Moving.y_median,...
            calibration_str.Fixed_binned_by_Moving.y_median - calibration_str.Fixed_binned_by_Moving.y_prctile(:, 3), ...
            calibration_str.Fixed_binned_by_Moving.y_prctile(:, 5) - calibration_str.Fixed_binned_by_Moving.y_median,...
            'LineStyle', 'none', 'Marker', 'x', 'LineWidth', 1.5);
    case {'Mean', 'mean'}
        errorbar(ax_hdl, calibration_str.Fixed_binned_by_Moving.x_bin_val, ...
            calibration_str.Fixed_binned_by_Moving.y_mean ,...
            calibration_str.Fixed_binned_by_Moving.y_std, calibration_str.Fixed_binned_by_Moving.y_std,...
            'LineStyle', 'none', 'Marker', 'x', 'LineWidth', 1.5);
    otherwise 
        error('Unrecognized statistics name');
end
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.XScale = 'log';
ax_hdl.YScale = 'log';
ax_hdl.XLabel.String = sprintf('%s radius (\\mum)', strrep(calibration_str.Moving_image_group, '_', ' '));
ax_hdl.YLabel.String = sprintf('%s radius (\\mum)', strrep(calibration_str.Fixed_image_group , '_', ' '));
ax_hdl.XLim = [0, calibration_str.Hist_edge_range(2)];
ax_hdl.YLim = [0, calibration_str.Hist_edge_range(2)];
ax_hdl.XTick = [0.5 : 0.5 : 2, 3, 4 : 2 : calibration_str.Hist_edge_range(2)];
ax_hdl.YTick = [0.5 : 0.5 : 2, 3, 4 : 2 : calibration_str.Hist_edge_range(2)];

cbar_hdl = colorbar(ax_hdl);
ax_hdl.ColorScale = 'log';
cbar_hdl.Label.String = 'Number of data points';
cbar_hdl.Limits(1) = 1;
ax_hdl.FontSize = 14;
hold(ax_hdl, 'on');
plot(ax_hdl, calibration_str.Hist_edge, calibration_str.Hist_edge, 'k-.', 'LineWidth', 1);

end