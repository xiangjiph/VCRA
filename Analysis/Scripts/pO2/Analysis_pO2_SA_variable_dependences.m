% Run Analysis_scaling_pO2_dist_vs_vessel_length_stat first.
r_cap = pO2_SA.avg_cap_r_um;
krogh_coeff = 1/2;
%%
fit_info_cell = cell(5, 0);
fit_info_cell(:, end+1) = {'pO2_data.pO2_lm.dt_mean', 'pO2_data.pO2_lm.pO2_mean', 'Average d_{ulm} (\mum)', 'Average u_{lm}', 'dulm_2_ulm'};
fit_info_cell(:, end+1) = {'pO2_data.dt_lm.dt_mean', 'pO2_data.pO2_lm.pO2_mean', 'Average d_{lm} (\mum)', 'Average u_{lm}', 'dlm_2_ulm'};
fit_info_cell(:, end+1) = {'pO2_data.dt_lm.dt_mean', 'pO2_data.dt_lm.pO2_mean', 'Average d_{lm} (\mum)', 'Average u_{dlm}', 'dlm_2_udlm'};
num_fit = size(fit_info_cell, 2);
%% Try to fit Krogh to local pO2 minimum vs its corresponding dt
scaled_krogh_str = struct;
for iter_fit = 1 : num_fit
    tmp_x_field = fit_info_cell{1, iter_fit};
    tmp_y_field = fit_info_cell{2, iter_fit};
    tmp_x_label = fit_info_cell{3, iter_fit};
    tmp_y_label = fit_info_cell{4, iter_fit};
    tmp_folder_name = fit_info_cell{5, iter_fit};
    
    tmp_x = cellfun(@(x) fun_getfield(x, tmp_x_field), wb_data_cell, 'UniformOutput', false);
    tmp_x = cat(1, tmp_x{:});
    tmp_y = cellfun(@(x) fun_getfield(x, tmp_y_field), wb_data_cell, 'UniformOutput', false);
    tmp_y = cat(1, tmp_y{:});
    tmp_fit_info_cell = cell(num_window_size, 1);
    for iter_wd_sz = 1 : num_window_size
        % Visualization for each window size
        tmp_wd_sz = lm_wz_um(iter_wd_sz);
        
        vis_x = tmp_x(:, iter_wd_sz);
        vis_y = tmp_y(:, iter_wd_sz);
        tmp_valid_Q = isfinite(vis_x) & isfinite(vis_y) & is_internal_cube_Q;
        tmp_num_data_before_selection = nnz(tmp_valid_Q);
        tmp_valid_Q = is_selected_cube_Q & tmp_valid_Q;
        vis_x = vis_x(tmp_valid_Q);
        vis_y = vis_y(tmp_valid_Q);
        fit_x = vis_x + r_cap;
        fit_y = - krogh_coeff * (fit_x .^ 2 .* (log(fit_x / r_cap) - 1/2) + r_cap^2/2);
        linear_fit_hdl = fitlm(fit_y, vis_y, 'Intercept', false);
        
        tmp_fit_info_str = struct;
        tmp_fit_info_str.window_size_um = tmp_wd_sz;
        tmp_fit_info_str.num_data_point = nnz(tmp_valid_Q);
        tmp_fit_info_str.fraction_of_valid_cube = tmp_fit_info_str.num_data_point / tmp_num_data_before_selection;
        tmp_fit_info_str.Estimate = linear_fit_hdl.Coefficients.Estimate.';
        tmp_fit_info_str.SE = linear_fit_hdl.Coefficients.SE.';
        tmp_fit_info_str.RSquaredAdj = linear_fit_hdl.Rsquared.Adjusted;
        tmp_fit_info_cell{iter_wd_sz} = tmp_fit_info_str;
        %% Output figures
        fig_hdl = figure('Visible', 'on');
        ax_hdl = axes(fig_hdl);
        histogram2(ax_hdl, vis_x, vis_y, 'DisplayStyle', 'tile');
        ax_hdl.XLim(1) = 0;
        ax_hdl.YLim(2) = 0;
        hold(ax_hdl, 'on');
        plot_x = 0 : ax_hdl.XLim(2);
        plot_fit_y = - krogh_coeff * ((plot_x + r_cap) .^ 2 .* (log((plot_x + r_cap) / r_cap) - 1/2) + r_cap^2/2);
        line_hdl = plot(ax_hdl, plot_x, plot_fit_y .* linear_fit_hdl.Coefficients.Estimate(1), ...
            'LineWidth', 2);
        ax_hdl.ColorScale = 'log';
        cbar_hdl = colorbar(ax_hdl);
        cbar_hdl.Limits(1) = 1;
        cbar_hdl.Label.String = 'Number of cubes';
        ax_hdl.FontSize = 14;
        ax_hdl.XLabel.String = tmp_x_label;
        ax_hdl.YLabel.String = tmp_y_label;
        ax_hdl.Title.String = sprintf('Window size %d \\mum', tmp_wd_sz);
        
        leg_str = sprintf('Scaled Krogh model\n\\lambda: %.4f \\pm %.4f \nR^2-Adjusted: %.3f\nData size: %d (%.3f)', ...
            linear_fit_hdl.Coefficients.Estimate(1), linear_fit_hdl.Coefficients.SE(1), ...
            linear_fit_hdl.Rsquared.Adjusted, ...
            tmp_fit_info_str.num_data_point, tmp_fit_info_str.fraction_of_valid_cube);
        leg_hdl = legend(ax_hdl, line_hdl, leg_str, 'Location', 'southwest');
        
        tmp_fp = fullfile(save_im_folder, tmp_folder_name, sprintf('%s_%s_%s_corr_coeff_min_c2v_vr_%.2f_wd_sz_%d_um.png', ...
            dataset_name, merge_stack_name, tmp_folder_name, ...
            min_cap2vsl_vol_fraction, tmp_wd_sz));
        %%
        fun_print_image_in_several_formats(fig_hdl, tmp_fp);
        delete(fig_hdl);
    end
    scaled_krogh_str.(tmp_folder_name) = cat(1, tmp_fit_info_cell{:});
end
%%
scaled_krogh_str.r_cap = r_cap;
scaled_krogh_str.krogh_coeff = krogh_coeff;
pO2_SA.scaled_krogh_coeff = scaled_krogh_str;
DataManager.write_data(pO2_SA.filepath, pO2_SA);
% The remaining questions are: why does the local extrema work? Because the
% boundary condition was effectively matched? Then why does the the attempt
% to fit the pO2 at other position with a scalar correction factor to the
% Krogh model did not work?
%% Bin pO2_vs_dt curve according to pO2 lm dt v
% Evaluate the pO2 value at different dt and check whether the dependence
% on pO2_lm_dt_v is parabolic
% Conclusion: Simple scaling of the Krogh model can not fit the curve well
%
% Get u(d) interpolation
pO2_in_dt_med = cellfun(@(x) x.pO2_data.pO2_stat_in_dt_bin.y_mean, wb_data_cell, 'UniformOutput', false);
pO2_in_dt_med = cat(1, pO2_in_dt_med{:});
pO2_dt_bin_val = cellfun(@(x) x.pO2_data.pO2_stat_in_dt_bin.x_bin_val, wb_data_cell, 'UniformOutput', false);
pO2_dt_bin_val = cat(1, pO2_dt_bin_val{:});
num_cube = numel(pO2_dt_bin_val);
pO2_vs_dt_itp = cell(num_cube, 1);
for iter_cube = 1 : num_cube
    tmp_x = pO2_dt_bin_val{iter_cube};
    tmp_y = pO2_in_dt_med{iter_cube};
    if ~isempty(tmp_x)
        pO2_vs_dt_itp{iter_cube} = griddedInterpolant(tmp_x, tmp_y, 'linear');
    end
end
%% Get binned cube list idx
iter_wd_sz = 11;
tmp_wd_sz = lm_wz_um(iter_wd_sz);

bin_data = cellfun(@(x) fun_getfield(x, 'pO2_data.pO2_lm.dt_mean'), wb_data_cell, 'UniformOutput', false);
bin_data = cat(1, bin_data{:});
bin_data = bin_data(:, iter_wd_sz);


bin_data_edge = 20 : 1 : 40;
bin_data_val = movmean(bin_data_edge, 2, 'Endpoints', 'discard');
bin_data_ind_cell = fun_bin_data_to_idx_list_by_edges(bin_data, bin_data_edge);
%% Try to fit the u(d) by parabolic 
% vis_bin_idx = 1 : 1 : (numel(bin_data_edge) - 1);
% vis_evaluate_dist_list = 5 : 5 : 35;
% num_vis_evaluate_dist = numel(vis_evaluate_dist_list);
% fig_hdl = figure;
% ax_hdl = axes(fig_hdl);
% [fit_info.k, fit_info.b, fit_info.RSquaredAdjusted]= deal(nan(num_vis_evaluate_dist, 1));
% line_hdl_array = [];
% for iter_evl = 1 : num_vis_evaluate_dist
%     vis_evaluate_dist = vis_evaluate_dist_list(iter_evl);
%     tmp_vis_bin_idx = vis_bin_idx;
%     tmp_vis_bin_idx = tmp_vis_bin_idx(bin_data_val >= vis_evaluate_dist);
%     vis_bin_val = bin_data_val(tmp_vis_bin_idx);
%     num_vis_bin = numel(tmp_vis_bin_idx);
%     vis_u0_cell = cell(num_vis_bin, 1);
%     for iter_cell = 1 : numel(tmp_vis_bin_idx)
%         tmp_idx = bin_data_ind_cell{tmp_vis_bin_idx(iter_cell)};
%         tmp_idx = tmp_idx(is_selected_cube_Q(tmp_idx));
%         tmp_u0 = cellfun(@(x) x(vis_evaluate_dist), ...
%             pO2_vs_dt_itp(tmp_idx));
%         tmp_u0 = tmp_u0(isfinite(tmp_u0));
%         vis_u0_cell{iter_cell} = tmp_u0;
%     end
%     vis_u0_vec = cat(1, vis_u0_cell{:});
%     vis_u0_class = repelem((1 : num_vis_bin)', cellfun(@numel, vis_u0_cell), 1);
%     vis_u0_prctile = cellfun(@(x) prctile(x, [25, 50, 75]), vis_u0_cell, 'UniformOutput', false);
%     vis_u0_prctile = cat(1, vis_u0_prctile{:});
%     
%     [ax_hdl, tmp_line_hdl, tmp_patch_hdl] = fun_vis_confidence_interval_shaded(vis_bin_val, ...
%         vis_u0_prctile(:, 2), vis_u0_prctile(:, 1), vis_u0_prctile(:, 3), ax_hdl);
%     line_hdl_array(end+1) = tmp_line_hdl;
%     hold(ax_hdl, 'on');
%     %% Try to fix a parabolic
%     fit_x = (vis_bin_val + r_cap) .^ 2;
%     fit_y = vis_u0_prctile(:, 2);
%     fit_hdl = fitlm(fit_x, fit_y, 'Intercept', true);
%     
%     fit_info.b(iter_evl) = fit_hdl.Coefficients.Estimate(1);
%     fit_info.k(iter_evl) = fit_hdl.Coefficients.Estimate(2);
%     fit_info.RSquaredAdjusted(iter_evl) = fit_hdl.Rsquared.Adjusted;
% end
% ax_hdl.XLabel.Interpreter = 'latex';
% ax_hdl.XLabel.String = 'Mean $d_{ulm}$ ($\mu m$)';
% ax_hdl.YLabel.Interpreter = 'latex';
% ax_hdl.YLabel.String = 'Median $u(d)$';
% leg_hdl = legend(ax_hdl, line_hdl_array, arrayfun(@(x) num2str(x, '%d'), vis_evaluate_dist_list, 'UniformOutput', false), ...
%     'Location', 'southwest');
% leg_hdl.Title.Interpreter = 'latex';
% leg_hdl.Title.String = '$d\;(\mu m)$';
% ax_hdl.FontSize = 14;

% tmp_fp = fullfile(save_im_folder, sprintf('%s_%s_u_vs_pO2_lm_dt_diff_d_cap2vsl_vol_gt_%03d.png', ...
%     dataset_name, merge_stack_name, round(min_cap2vsl_vol_fraction * 100)));
% fun_print_image_in_several_formats(fig_hdl, tmp_fp);
% delete(fig_hdl);
%% Try to fit Krogh
vis_bin_idx = 1 : 1 : (numel(bin_data_edge) - 1);
vis_r_list = 5 : 5 : 35;
num_vis_evaluate_dist = numel(vis_r_list);
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.2;
ax_hdl = axes(fig_hdl);
[fit_info.k, fit_info.b, fit_info.RSquaredAdjusted]= deal(nan(num_vis_evaluate_dist, 1));
line_hdl_array = [];
fit_line_hdl_array = [];
for iter_evl = 1 : num_vis_evaluate_dist
    vis_r = vis_r_list(iter_evl);
    tmp_vis_bin_idx = vis_bin_idx;
    % Determine the cube with average dulm greater than r being evaluated
    tmp_vis_bin_idx = tmp_vis_bin_idx(bin_data_val >= vis_r);
    vis_bin_val = bin_data_val(tmp_vis_bin_idx);
    num_vis_bin = numel(tmp_vis_bin_idx);
    % Evaluate the median u(d) for all the 240 cubes in each bin
    vis_u0_cell = cell(num_vis_bin, 1);
    for iter_cell = 1 : numel(tmp_vis_bin_idx)
        tmp_idx = bin_data_ind_cell{tmp_vis_bin_idx(iter_cell)};
        tmp_idx = tmp_idx(is_selected_cube_Q(tmp_idx));
        tmp_u0 = cellfun(@(x) x(vis_r), ...
            pO2_vs_dt_itp(tmp_idx));
        tmp_u0 = tmp_u0(isfinite(tmp_u0));
        vis_u0_cell{iter_cell} = tmp_u0;
    end
    vis_u0_prctile = cellfun(@(x) prctile(x, [25, 50, 75]), vis_u0_cell, 'UniformOutput', false);
    vis_u0_prctile = cat(1, vis_u0_prctile{:});
    
    line_med_hfl = plot(ax_hdl, vis_bin_val, vis_u0_prctile(:, 2), 'LineWidth', 2);
    line_hdl_array(end+1) = line_med_hfl;
    hold(ax_hdl, 'on');
    %% Try to fit the Krogh model
    fit_x = - krogh_coeff * ((vis_bin_val + r_cap) .^ 2 .* log((vis_r + r_cap)/ r_cap) ...
        - ((vis_r + r_cap)^2 - r_cap^2)/2);
    fit_y = vis_u0_prctile(:, 2);
    fit_hdl = fitlm(fit_x, fit_y, 'Intercept', false);
    line_hdl = plot(ax_hdl, vis_bin_val, fit_x .* fit_hdl.Coefficients.Estimate(1), ...
        'LineStyle', ':', 'LineWidth', 2, 'Color', line_med_hfl.Color);
    fit_line_hdl_array(end+1) = line_hdl;
    fit_info.b(iter_evl) = fit_hdl.Coefficients.Estimate(1);
    fit_info.RSquaredAdjusted(iter_evl) = fit_hdl.Rsquared.Adjusted;
end
ax_hdl.XLabel.Interpreter = 'latex';
ax_hdl.XLabel.String = 'Mean $d_{ulm}$ ($\mu m$)';
ax_hdl.YLabel.Interpreter = 'latex';
ax_hdl.YLabel.String = 'Median $u(d)$';
ax_hdl.FontSize = 14;
leg_hdl = legend(ax_hdl, fit_line_hdl_array, ...
    arrayfun(@(x1, x2, x3) sprintf('d = %02d \\mum, \\lambda = %.2f, R^2 = %.2f', x1, x2, x3),...
    vis_r_list', fit_info.b, fit_info.RSquaredAdjusted, 'UniformOutput', false), ...
    'Location', 'southwest');
leg_hdl.Title.String = 'Fitting';
leg_hdl.FontSize = 10;

% tmp_fp = fullfile(save_im_folder, sprintf('%s_%s_u_vs_ulm_dt_scale_Krogh_cap2vsl_vol_gt_%03d.png', ...
%     dataset_name, merge_stack_name, round(min_cap2vsl_vol_fraction * 100)));
% fun_print_image_in_several_formats(fig_hdl, tmp_fp);
% delete(fig_hdl);