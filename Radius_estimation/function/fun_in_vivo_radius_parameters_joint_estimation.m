function est_para_str = fun_in_vivo_radius_parameters_joint_estimation(psf_est_int, tile_image, skel_sub, ...
    vsl_r_est_0, ori_vec, voxel_size, radius_estimation_method)
%% 
% Input: 
%   skel_sub: N-by-3 integer array. The default unit of skel_sub is
%   micron.The skeleton voxel subscripts would not be rescaled inside this
%   function. 
%%
verbose_Q = false;
visualization_profile_Q = false;
%% Initial estimation of PSF
num_round_est = 3;
num_r_est_round = 6;

est_para_str = [];
best_fit_corr = -1;

med_skel_sub = skel_sub(ceil(size(skel_sub, 1)/2), :);

tile_im_size = size(tile_image);
min_dist_to_boundary_z = ceil(vsl_r_est_0 * 2);
% Not enough information for fitting the z-profile
if med_skel_sub(3) < min_dist_to_boundary_z || ...
        med_skel_sub(3) + min_dist_to_boundary_z > tile_im_size(3)
    return;
end
%%
ori_vec_z = abs(ori_vec(3));
for iter_est = 1 : num_round_est 
    if isnan(vsl_r_est_0)
        est_para_str = [];
        return;
    end
    %% Get radial intensity profile
    % Interpolate local on-plane intensity profile
    int_profile_str_xy = fun_radius_estimation_get_xy_int_profile(tile_image, ...
        med_skel_sub, vsl_r_est_0, ori_vec, voxel_size);
    
    int_profile_str_z = fun_radius_estimation_get_z_int_profile(tile_image, ...
        med_skel_sub, vsl_r_est_0, ori_vec(3), int_profile_str_xy.min_int, ...
        voxel_size);
    %% Generate intensity profiles for each PSF
    % Find the closest simulated vessel parameters
    est_para_str = fun_in_vivo_radius_estimation_get_best_fit_parameters(psf_est_int, vsl_r_est_0, ...
        ori_vec_z, int_profile_str_xy, int_profile_str_z);
    assert(isfinite(est_para_str.corr_max), 'Angle-weighted correlation coefficient is not finite'); 
    %%
    if est_para_str.corr_max > best_fit_corr
        best_fit_corr = est_para_str.corr_max;   
        best_fit_psf_str = est_para_str;       
        % Update radius estimation
        % Use multiple algorithm to estimate the radius
        % Distance transform
        n_edge_int_itp = psf_est_int.edge_int_str{est_para_str.best_fit_rbc_int_n_ind, ...
            est_para_str.best_fit_blt_ind}.n_min_edge_int_interpolation.n_min_edge_int;
        vsl_r_est_dt = fun_in_vivo_radius_estimation_by_DT(tile_image, skel_sub, vsl_r_est_0, ori_vec_z, ...
            n_edge_int_itp, voxel_size, num_r_est_round);
        
        % Profile with calculated threshold
        n_edge_int_itp = psf_est_int.edge_int_str{est_para_str.best_fit_rbc_int_n_ind, ...
            est_para_str.best_fit_blt_ind}.max_n_min_edge_int_interpolation.n_min_edge_int;
        vsl_r_est_prof_c = fun_radius_estimation_by_xy_int_profile(int_profile_str_xy, vsl_r_est_0, ...
            ori_vec_z, n_edge_int_itp, num_r_est_round);
        
        % Profile with normalized edge intensity 0.5
        vsl_r_est_prof_half = fun_radius_estimation_by_xy_int_profile(int_profile_str_xy, vsl_r_est_0, ...
            ori_vec_z, 0.5, num_r_est_round);
        
        switch radius_estimation_method
            case 'DT'
                vsl_r_est_0 = vsl_r_est_dt;
            case 'Prof'
                vsl_r_est_0 = vsl_r_est_prof_c;
            case 'Prof_half'
                vsl_r_est_0 = vsl_r_est_prof_half;
        end
        if isnan(vsl_r_est_0)
            break;
        end
    else
%         fprintf('New estimation does not increase the fitting accuracy. Return\n');
%         if ori_vec_z < 1/2 && best_fit_corr > 0.99 && iter_est >= 2
%             fprintf('Visualization\n');
%         end
        break;
    end
    %% Visualization
    if visualization_profile_Q
        vis_tick_step = 10;
        fig_hdl = figure;
        fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 2];
        ax_hdl_1 = subplot(2,3,1);
        imagesc(ax_hdl_1, best_fit_psf_str.corr_mat);
        colormap(ax_hdl_1, 'jet');
        cbar_hdl = colorbar(ax_hdl_1);
        cbar_hdl.Label.String = 'Correlation';
        ax_hdl_1.XTick = 1 : vis_tick_step : numel(psf_est_int.vsl_bc_lyr_thkn_um_list);
        ax_hdl_1.XTickLabel = arrayfun(@(x) num2str(x, '%.2f'), psf_est_int.vsl_bc_lyr_thkn_um_list(1:vis_tick_step:end), ...
            'UniformOutput', false);
        ax_hdl_1.YTick = 1 : vis_tick_step : numel(psf_est_int.rbc_int_n_list);
        ax_hdl_1.YTickLabel = arrayfun(@(x) num2str(x, '%.2f'), psf_est_int.rbc_int_n_list(1:vis_tick_step:end), ...
            'UniformOutput', false);
        ax_hdl_1.XLabel.String = 'BLT (\mum)';
        ax_hdl_1.YLabel.String = 'RBC Intensity';
        ax_hdl_1.Title.String = sprintf('Avg. Corr. (max = %.5f)', best_fit_psf_str.corr_max);
        
        ax_hdl_2 = subplot(2,3,2);
        imagesc(ax_hdl_2, best_fit_psf_str.xy_corr_mat);
        colormap(ax_hdl_2, 'jet');
        cbar_hdl = colorbar(ax_hdl_2);
        cbar_hdl.Label.String = 'Correlation';
        ax_hdl_2.XTick = 1 : vis_tick_step : numel(psf_est_int.vsl_bc_lyr_thkn_um_list);
        ax_hdl_2.XTickLabel = arrayfun(@(x) num2str(x, '%.2f'), psf_est_int.vsl_bc_lyr_thkn_um_list(1:vis_tick_step:end), ...
            'UniformOutput', false);
        ax_hdl_2.YTick = 1 : vis_tick_step : numel(psf_est_int.rbc_int_n_list);
        ax_hdl_2.YTickLabel = arrayfun(@(x) num2str(x, '%.2f'), psf_est_int.rbc_int_n_list(1:vis_tick_step:end), ...
            'UniformOutput', false);
        ax_hdl_2.XLabel.String = 'BLT (\mum)';
        ax_hdl_2.YLabel.String = 'RBC Intensity';
        ax_hdl_2.Title.String = 'Corr. XY';
        
        ax_hdl_3 = subplot(2,3,3);
        imagesc(ax_hdl_3, best_fit_psf_str.z_corr_mat);
        colormap(ax_hdl_3, 'jet');
        cbar_hdl = colorbar(ax_hdl_3);
        cbar_hdl.Label.String = 'Correlation';
        ax_hdl_3.XTick = 1 : vis_tick_step : numel(psf_est_int.vsl_bc_lyr_thkn_um_list);
        ax_hdl_3.XTickLabel = arrayfun(@(x) num2str(x, '%.2f'), psf_est_int.vsl_bc_lyr_thkn_um_list(1:vis_tick_step:end), ...
            'UniformOutput', false);
        ax_hdl_3.YTick = 1 : vis_tick_step : numel(psf_est_int.rbc_int_n_list);
        ax_hdl_3.YTickLabel = arrayfun(@(x) num2str(x, '%.2f'), psf_est_int.rbc_int_n_list(1:vis_tick_step:end), ...
            'UniformOutput', false);
        ax_hdl_3.XLabel.String = 'BLT (\mum)';
        ax_hdl_3.YLabel.String = 'RBC Intensity';
        ax_hdl_3.Title.String = sprintf('Corr. Z (v_z = %.2f)', ori_vec_z);
        
        ax_hdl_4 = subplot(2,3,4);
        plot(ax_hdl_4, int_profile_str_xy.r, int_profile_str_xy.int_val_n_smooth);
        hold(ax_hdl_4, 'on');
        plot(ax_hdl_4, best_fit_psf_str.xy_max_corr_int_r + best_fit_psf_str.xy_max_corr_shift_um, best_fit_psf_str.xy_max_corr_int_profile);
        ax_hdl_4.XLabel.String = 'r (\mum)';
        ax_hdl_4.YLabel.String = 'Normalized intensity';
        legend('Measurement', 'Best fit');
        ax_hdl_4.Title.String = sprintf('Lateral intensity profile\nr_{est} = %.1f \\mum, RBC Int_n %.1f, BLT %.1f \\mum', ...
            vsl_r_est_0, best_fit_psf_str.best_fit_parameter);
        
        ax_hdl_5 = subplot(2,3,5);
        plot(ax_hdl_5, int_profile_str_z.r, int_profile_str_z.int_val_n_smooth);
        hold(ax_hdl_5, 'on');
        plot(ax_hdl_5, best_fit_psf_str.z_max_corr_int_r + best_fit_psf_str.z_max_corr_shift_um, best_fit_psf_str.z_max_corr_int_profile);
        ax_hdl_5.XLabel.String = 'r (\mum)';
        ax_hdl_5.YLabel.String = 'Normalized intensity';
        legend('Measurement', 'Best fit');
        ax_hdl_5.Title.String = sprintf('Axial intensity profile\nr_{est} = %.1f \\mum', vsl_r_est_0);
        %% Save visualization
        dataset_name = 'Vessel_radius_calibration';
        stack = 'In_vivo_image_simulation';
        DataManager = FileManager;
        im_save_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), ...
            'Para_joint_est');
        fig_fp = fullfile(im_save_folder, sprintf('In_vivo_im_para_jt_est_vsl_r_est_%d_nm_z_comp_%d.png', ...
            round(1000 * vsl_r_est_0), round(ori_vec_z * 100)));
        fun_print_image_in_several_formats(fig_hdl, fig_fp);
    end
end
est_para_str = best_fit_psf_str;
est_para_str.radius_est_method = {'DT', 'Profile_computed_int', 'Profile_half_int'};
est_para_str.vsl_r_est = [vsl_r_est_dt, vsl_r_est_prof_c, vsl_r_est_prof_half];
est_para_str.vsl_ori_vec = ori_vec;
end