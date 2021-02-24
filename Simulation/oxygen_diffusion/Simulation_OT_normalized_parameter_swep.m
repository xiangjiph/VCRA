%% Load capillary network
set_env;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
load_skel_ver = '240_cube_re';
grid_info = DataManager.load_grid(dataset_name, stack, '240_cube');
mask_version = '240_cube_recon';
vis_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'Tissue_pO2_simulation');
%%
inhomogeneous_term = 1;
krogh_coeff = 1/2;
r_cap = 2;
cube_label_list = [24119, 24148, 24147, 24174, 24204, 24262];
num_cubes = numel(cube_label_list);
[pO2_result, pO2_analyzed_result] = deal(cell(num_cubes, 1));
[cap_dt_max_med, local_dt_max_med, local_oxy_min_dt_med] = deal(nan(num_cubes, 1));
downsample_rate = 2;
visualization_Q = false;
save_figure_Q = true;
parfor iter_cube = 1 : num_cubes
    tmp_cube_label = cube_label_list(iter_cube);
    % tmp_cube_label = 31998;
    tmp_cube_grid_sub = grid_info.bbox_grid_sub_list(tmp_cube_label, :);
    tmp_recon_str = DataManager.load_block_mask(dataset_name, stack, ...
        mask_version, tmp_cube_grid_sub(1), tmp_cube_grid_sub(2), tmp_cube_grid_sub(3));
    tmp_recon = fun_reconstruct_block_mask(tmp_recon_str);
    cap_dt_max_med(iter_cube) = median(tmp_recon_str.link.features.nearest_tissue_dt_max, 'omitnan');
    % volumeViewer(vis_recon);
    tmp_target_size = round(size(tmp_recon) / downsample_rate);
    tmp_recon_rz = imresize3(uint8(tmp_recon), tmp_target_size) > 0;
    mask_rz_dt = bwdist(tmp_recon_rz) .* downsample_rate;
    tmp_inhomgeneous_term = inhomogeneous_term * downsample_rate.^ 2;
    %% Use Krogh model to estimate the initial condition
    if isfinite(cap_dt_max_med(iter_cube))
        local_median_dt_max = cap_dt_max_med(iter_cube);
        fit_x_fun = @(x) (- krogh_coeff * (local_median_dt_max .^ 2 .* log((x + r_cap) ./ r_cap) - ...
            ((x + r_cap).^ 2 - r_cap^2)/2));
        tmp_ini_pO2_n = 0.5 * fit_x_fun(mask_rz_dt);
    else
        tmp_ini_pO2_n = -300;
    end
    pO2_n_result_ds2 = fun_simulation_OT_solve_ct_diff_itr_n(tmp_recon_rz, tmp_inhomgeneous_term, tmp_ini_pO2_n, true);  
    pO2_n_result_ds2.vessel_mask_dt = mask_rz_dt;
    
    [oxy_result_n_str, tmp_x_data, tmp_y_data] = fun_simulation_OT_analyze_result_n(...
        pO2_n_result_ds2.vessel_mask_dt, pO2_n_result_ds2.pO2_array, 16, 10);
    pO2_result{iter_cube} = pO2_n_result_ds2;
    pO2_analyzed_result{iter_cube} = oxy_result_n_str;
    local_dt_max_med(iter_cube) = median(oxy_result_n_str.dt_lm.dt_v, 'omitnan');
    local_oxy_min_dt_med(iter_cube) = median(oxy_result_n_str.pO2_lm.dt_v, 'omitnan');
    %% Use krogh model to estimate initial pO2 distribution    
    fit_x = fit_x_fun(tmp_x_data);
    lin_fit_hdl = fitlm(fit_x, tmp_y_data, 'Intercept', false);
    correction_factor = lin_fit_hdl.Coefficients.Estimate(1);
    %% Visualization
    if visualization_Q
        fig_hdl = figure;
        ax_hdl = axes(fig_hdl);
        fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
        histogram2(ax_hdl, tmp_x_data, tmp_y_data, 'DisplayStyle', 'tile');
        hold(ax_hdl, 'on');
        % med_hdl = plot(ax_hdl, tmp_y_stat.x_bin_val, tmp_y_stat.y_median, 'LineWidth', 2);
        mean_hdl = plot(ax_hdl, oxy_result_n_str.pO2_stat_in_bin.x_bin_val, oxy_result_n_str.pO2_stat_in_bin.y_mean, 'LineWidth', 4);
        ptl005_hdl = plot(ax_hdl, oxy_result_n_str.pO2_stat_in_bin.x_bin_val, oxy_result_n_str.pO2_stat_in_bin.y_prctile(:, 1), 'LineWidth', 4, ...
            'LineStyle', '-.');
        ptl095_hdl = plot(ax_hdl, oxy_result_n_str.pO2_stat_in_bin.x_bin_val, oxy_result_n_str.pO2_stat_in_bin.y_prctile(:, end), 'LineWidth', 4, ...
            'LineStyle', '-.');
        % correction_factor = 1;
        plot_krogh_x = 0 : min(local_median_dt_max, max(tmp_x_data));
        oxygen_krogh = correction_factor * fit_x_fun(plot_krogh_x);
        plt_krogh = plot(ax_hdl, plot_krogh_x, oxygen_krogh, 'LineWidth', 4, 'LineStyle', '-');
        
        oxygen_krogh_0 = fit_x_fun(plot_krogh_x);
        plt_krogh_0 = plot(ax_hdl, plot_krogh_x, oxygen_krogh_0, 'LineWidth', 4, 'LineStyle', '-');
        
        dot_dt_lm_hdl = scatter(ax_hdl, oxy_result_n_str.dt_lm.dt_v, oxy_result_n_str.dt_lm.pO2_v, 'y+');
        dot_oxy_lm_hdl = scatter(ax_hdl, oxy_result_n_str.pO2_lm.dt_v, oxy_result_n_str.pO2_lm.pO2_v, 'r*');
        %
        leg_hdl_array = [mean_hdl, ptl005_hdl, ptl095_hdl, plt_krogh, plt_krogh_0, dot_dt_lm_hdl, dot_oxy_lm_hdl];
        leg_hdl_string = {'Average', '5% percentile', '95% percentile', ...
            sprintf('Fitting Krogh model:\nCorr. factor: %.2f\nR^2-Adjusted: %.2f', correction_factor, ...
            lin_fit_hdl.Rsquared.Adjusted), 'Krogh model'...
            'Local d maximum', 'Local u_0 minimum'};
        legend(ax_hdl, leg_hdl_array, leg_hdl_string);
        ax_hdl.XLabel.String = 'd';
        ax_hdl.YLabel.String = 'u_0(d)';
        ax_hdl.FontSize = 14;
        if save_figure_Q
            fig_fp = fullfile(vis_folder, sprintf('%s_%s_compare_numerical_to_Krogh_model_240_cube_label_%d.png', ...
                dataset_name, stack, cube_label_list(iter_cube)));
            fun_print_image_in_several_formats(fig_hdl, fig_fp);
        end        
    end    
end
%% Compare the numerical curve in different regions
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
for iter_cube = 1 : num_cubes
    plot(ax_hdl, pO2_analyzed_result{iter_cube}.pO2_stat_in_bin.x_bin_val, ...
        pO2_analyzed_result{iter_cube}.pO2_stat_in_bin.y_mean, 'LineWidth', 2);
    hold(ax_hdl, 'on');
end
