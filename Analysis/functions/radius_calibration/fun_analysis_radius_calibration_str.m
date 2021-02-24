function tmp_str = fun_analysis_radius_calibration_str(fixed_r, moving_r, ...
    vis_stat_name, hist_edge_range, num_hist_bin, flat_r_max, equal_r_min, spline_skip_step, fit_min_num_data)

tmp_str = struct;
tmp_str.Fixed_r_um = fixed_r;
tmp_str.Moving_r_um = moving_r;
tmp_str.Hist_edge_range = hist_edge_range;
tmp_str.Hist_edge = 10 .^ [linspace(log10(hist_edge_range(1)), ...
    log10(hist_edge_range(2)), num_hist_bin + 1)];
tmp_str.Fixed_binned_by_Moving = fun_analysis_get_y_stat_in_x_bin(tmp_str.Moving_r_um , tmp_str.Fixed_r_um, tmp_str.Hist_edge);
tmp_str.Num_segments = numel(tmp_str.Fixed_r_um);
% Interpolation
itp_x = tmp_str.Fixed_binned_by_Moving.x_bin_val;
itp_y = tmp_str.Fixed_binned_by_Moving.(sprintf('y_%s', vis_stat_name));
itp_y_flat_value = mean(itp_y(find(itp_x >= flat_r_max, 2, 'first')));
itp_y(itp_x < flat_r_max) = itp_y_flat_value;
kept_Q = (itp_x < equal_r_min);
kept_Q(2 : spline_skip_step : end) = false;
equla_r_point_start = find(kept_Q, 2, 'last');
equla_r_point_start = abs(diff(itp_x(equla_r_point_start)));
itp_y = itp_y(kept_Q);
itp_x = itp_x(kept_Q);
equal_r_points = [(itp_x(end) + equla_r_point_start) : 1 : 100]';
itp_y = cat(1, itp_y, equal_r_points);
itp_x = cat(1, itp_x, equal_r_points);
% Zero the negative 
% pad_x = [-4 : 1 : 0, hist_edge_range(1) - (5 : -1 : 1) * (10 * eps)]';
% pad_y = [zeros(1, 6), itp_y_flat_value * ones(1, 4)]';
% itp_x = cat(1, pad_x, itp_x);
% itp_y = cat(1, pad_y, itp_y);

% Add asymptotic points
tmp_str.spline_itp = griddedInterpolant(itp_x, itp_y, 'spline', 'linear');
% Guess an analytical function
tmp_fit_x = 1 ./ (1 + tmp_str.Fixed_binned_by_Moving.x_bin_val);
tmp_fit_y = tmp_str.Fixed_binned_by_Moving.(sprintf('y_%s', vis_stat_name)) - tmp_str.Fixed_binned_by_Moving.x_bin_val;
fit_selected_Q = tmp_str.Fixed_binned_by_Moving.y_count > fit_min_num_data;
tmp_fit_x = tmp_fit_x(fit_selected_Q);
tmp_fit_y = tmp_fit_y(fit_selected_Q);
tmp_fit = fitlm(tmp_fit_x, tmp_fit_y, 'Intercept', false);
tmp_str.formula_r0 = tmp_fit.Coefficients.Estimate(1);
tmp_str.formula_R2 = tmp_fit.Rsquared.Adjusted;
tmp_str.formula = @(x) (x + tmp_str.formula_r0 ./ (1 + x));
end