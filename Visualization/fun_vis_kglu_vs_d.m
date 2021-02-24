function [fig_hdl, varargout]= fun_vis_kglu_vs_d(data_table, vis_wz_idx, x_stat_name, ...
    x_label, stack_list, cap_r_um, scaled_Krogh_coeff_str)
%% Parameters
glucoseOxygenIndex = 5.65;
alphaO2 = 9.8214e-4 / 760;
% Oxygen diffusion coefficient
% 2.78e-9 m^2/s in water at 40.2 C
% 2.52e-9 m^2/s in water at 35.1 C
diffCoeffO2 = 1.9e-9; % m^2/s in rodent brain (Clark et al 1978)
% Oxygen solubility
% 1.3e-3 mol(L * atm) in water;
% 9.8214e-4 mol/(L * atm) in small rodent brain (in vivo) Clark et al 1978
% alphaO2 = 1.3e-3 / 760; % mol/(L * mmHg) in water
%%
tmp_data_x = data_table.(sprintf('%s_%d_median', x_stat_name, vis_wz_idx));
tmp_data_x_25_add = data_table.(sprintf('%s_%d_ptl025', x_stat_name, vis_wz_idx)) - tmp_data_x;
tmp_data_x_75_add = data_table.(sprintf('%s_%d_ptl075', x_stat_name, vis_wz_idx)) - tmp_data_x;

num_stack = numel(stack_list);
tmp_data_y_mean = data_table.k_glu_mean_M_per_sec * 1e6;
tmp_data_y_std = data_table.k_glu_std_M_per_sec * 1e6;
tmp_fit_y = repmat(tmp_data_y_mean, 1, num_stack);
tmp_fit_y_std = repmat(tmp_data_y_std, 1, num_stack);

scale_krogh_coeff = scaled_Krogh_coeff_str(vis_wz_idx).Estimate;
vis_wz_um = scaled_Krogh_coeff_str(vis_wz_idx).window_size_um;

plot_residual_Q = true;
%% Add selection here
if ismember('included_for_fitting_Q', data_table.Properties.VariableNames)
    selected_Q = data_table.included_for_fitting_Q;
else
    selected_Q = true(size(tmp_data_x));
end
tmp_fit_x = tmp_data_x(selected_Q);
tmp_fit_y = tmp_fit_y(selected_Q);
tmp_fit_y_std = tmp_fit_y_std(selected_Q);
%% Fitting
tmp_fit_fx_fun = @(x) (alphaO2 * diffCoeffO2 * 1e18) ./ (0.5 * glucoseOxygenIndex * scale_krogh_coeff *...
    ((x + cap_r_um) .^ 2 .* (log((x + cap_r_um) ./ cap_r_um) - 1/2) + cap_r_um ^ 2 / 2));
tmp_fit_fx = tmp_fit_fx_fun(tmp_fit_x);
linear_fit_hdl = fitlm(tmp_fit_fx, tmp_fit_y, 'Intercept', false);

avg_pO2_drop = linear_fit_hdl.Coefficients.Estimate(1);

tmp_fit_info_str = struct;
tmp_fit_info_str.window_size_um = vis_wz_um;
tmp_fit_info_str.num_region = numel(tmp_data_y_mean);
tmp_fit_info_str.Estimate = linear_fit_hdl.Coefficients.Estimate.';
tmp_fit_info_str.SE = linear_fit_hdl.Coefficients.SE.';
tmp_fit_info_str.RSquaredAdj = linear_fit_hdl.Rsquared.Adjusted;
if nargout > 1
    varargout{1} = tmp_fit_info_str;
end

best_fit_y = avg_pO2_drop .* tmp_fit_fx;
fit_residual = best_fit_y - tmp_fit_y;
fit_relative_residual = fit_residual ./ tmp_fit_y;
%% Visualization
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
if plot_residual_Q
    ax_hdl = subplot(4, 4, 1:12);
    ax_hdl_2 = subplot(4, 4, 13 : 16);
else
    ax_hdl = axes(fig_hdl);
end

for iter_stack = 1 : num_stack
    errorbar(ax_hdl, tmp_data_x(:, iter_stack), ...
        tmp_data_y_mean, tmp_data_y_std, tmp_data_y_std, tmp_data_x_25_add(:, iter_stack), tmp_data_x_75_add(:, iter_stack), 'LineStyle', 'none', 'LineWidth', 2);
    hold(ax_hdl, 'on');
end
tmp_plot_x = linspace(min(tmp_fit_x), max(tmp_fit_x), 20);
theo_plt_hdl = plot(ax_hdl, tmp_plot_x, avg_pO2_drop * tmp_fit_fx_fun(tmp_plot_x), ...
    'LineWidth', 3, 'Color', 'k');
box(ax_hdl, 'on');
ax_hdl.XLim(1) = 0;
ax_hdl.YLim(1) = 0;
ax_hdl.YLim(2) = 5 * ceil(ax_hdl.YLim(2)/5);
pO2_contour_val = [1, 5, 10, 20, 30, 40, 50, 60];
contour_plot_x_range = ax_hdl.XLim(1) : 0.5 : ax_hdl.XLim(2);
num_contour = numel(pO2_contour_val);
pO2_contour_cmap = jet(num_contour);
for iter_contour = 1 : num_contour
    tmp_contour_x = contour_plot_x_range;
    tmp_contour_y = pO2_contour_val(iter_contour) .* tmp_fit_fx_fun(contour_plot_x_range);
    plot(ax_hdl, tmp_contour_x, tmp_contour_y, ...
        'LineWidth', 1, 'Color', pO2_contour_cmap(iter_contour, :));
end
ax_hdl.YLabel.Interpreter = 'latex';
ax_hdl.YLabel.String = '$k_{glu} \;\frac{\mu mol}{s\cdot L}$';
legend(ax_hdl, {stack_list{:}, ...
    sprintf('Window size: %d \\mum\n\\DeltapO_2: %.1f \\pm %.1f mmHg\nR^2-Adjusted: %.3f\nNumber of regions: %d', ...
    vis_wz_um, avg_pO2_drop, linear_fit_hdl.Coefficients.SE(1), ...
    linear_fit_hdl.Rsquared.Adjusted, size(tmp_data_y_mean, 1))},...
    'Location', 'southwest', 'FontSize', 14);
ax_hdl.FontSize = 14;
ax_hdl.YLabel.FontSize = 18;

if plot_residual_Q
    errorbar(ax_hdl_2, tmp_fit_x, fit_relative_residual, tmp_fit_y_std ./ tmp_fit_y, 'LineStyle', 'none');
    grid(ax_hdl_2, 'on');
    max_abs_resdual = max(abs(fit_relative_residual));
    if max_abs_resdual <= 0.5
        ax_hdl_2.YTick = [-0.5 : 0.25 : 0.5];
        ax_hdl_2.YLim = [-0.5, 0.5];
    elseif max_abs_resdual <=1
        ax_hdl_2.YTick = [-1 : 0.5 : 1];
        ax_hdl_2.YLim = [-1, 1];
    end
    
    ax_hdl_2.XLim = ax_hdl.XLim;
    ax_hdl_2.YLabel.String = 'Relative error';
    if contains(x_label, '$')
        ax_hdl_2.XLabel.Interpreter = 'latex';
    end
    ax_hdl_2.XLabel.String = x_label;
    
    legend(ax_hdl_2, sprintf('Average absolute relative error: %.3f', mean(abs(fit_relative_residual))), 'Location', 'northwest');
    ax_hdl_2.FontSize = 14;
end
end