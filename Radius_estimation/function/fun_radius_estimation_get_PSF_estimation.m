function est_psf_str = fun_radius_estimation_get_PSF_estimation(psf_est_int, tile_image, skel_sub, ...
    vsl_r_est, ori_vec, voxel_size, psf_est_um, verbose_Q)
%% 
% Input: 
%   skel_sub: N-by-3 integer array. The default unit of skel_sub is
%   micron.The skeleton voxel subscripts would not be rescaled inside this
%   function. 
%%
if nargin < 7
    % The following PSF parameters are for Mouselight TP microscope
    psf_est_um = [0.8, 6.5];
    verbose_Q = false;
elseif nargin < 8
    verbose_Q = false;
end
visualization_profile_Q = false;
%% Initial estimation of PSF
initial_est_PSF_FWHM_xy = psf_est_um(1); % um
initial_est_PSF_FWHM_z = psf_est_um(2); % um
% The following PSF parameters are for DKLabe TP microscope
% initial_est_PSF_FWHM_xy = 0.6; % um
% initial_est_PSF_FWHM_z = 4.0; % um

[~, initial_est_PSF_FWHM_xy_ind] = min(abs(psf_est_int.psf_FWHM_xy_list - initial_est_PSF_FWHM_xy));
[~, initial_est_PSF_FWHM_z_ind] = min(abs(psf_est_int.psf_FWHM_z_list - initial_est_PSF_FWHM_z));
initial_est_PSF_list_ind = sub2ind(size(psf_est_int.psf_FWHM), initial_est_PSF_FWHM_xy_ind, initial_est_PSF_FWHM_z_ind);

est_psf_int_interpolation = psf_est_int.edge_int_str{initial_est_PSF_list_ind}.n_min_edge_int_interpolation.n_min_edge_int;
num_round_est = 3;
num_r_est_round = 6;

est_psf_str = [];
best_fit_corr = -1;

med_skel_sub = skel_sub(ceil(size(skel_sub, 1)/2), :);

tile_im_size = size(tile_image);
min_dist_to_boundary_z = ceil(vsl_r_est * 2);
% Not enough information for fitting the z-profile
if med_skel_sub(3) < min_dist_to_boundary_z || ...
        med_skel_sub(3) + min_dist_to_boundary_z > tile_im_size(3)
    return;
end
%% Determine the local intensity maximum
[local_im, local_im_bbox_mmxx] = crop_center_box(tile_image, med_skel_sub, [4,4,1]);
[local_int_max, local_int_max_ind] = max(local_im(:)); %#ok<ASGLU>
local_int_max_sub = fun_ind2sub(local_im_bbox_mmxx(4:6) - local_im_bbox_mmxx(1:3) + 1, ...
    local_int_max_ind);
tmp_int_max_sub = local_int_max_sub + local_im_bbox_mmxx(1:3) - 1;
%%
ori_vec_z = abs(ori_vec(3));
for iter_est = 1 : num_round_est 

    %% Estimate the radius
    [vsl_r_est, dt_max_sub] = fun_radius_estimation_by_adp_thrld_DT(tile_image, skel_sub, vsl_r_est, ori_vec, ...
        est_psf_int_interpolation, voxel_size, num_r_est_round);
    % The simulation is done for vessel of radius less than 5 um
    if isnan(vsl_r_est) || vsl_r_est > 5
        est_psf_str = [];
        return;
    end
    if verbose_Q
        fprintf('Estimated vessel radius is %f um, orientation vector z-component %f\n', vsl_r_est, ori_vec_z);
    end
    %% Get radial intensity profile
    % Interpolate local on-plane intensity profile
%     int_profile_str_xy = fun_radius_estimation_get_xy_int_profile(tile_image, ...
%         tmp_int_max_sub, vsl_r_est, ori_vec, voxel_size);
%     
%     int_profile_str_z = fun_radius_estimation_get_z_int_profile(tile_image, ...
%         tmp_int_max_sub, vsl_r_est, ori_vec(3), int_profile_str_xy.min_int, ...
%         voxel_size);
    int_profile_str_xy = fun_radius_estimation_get_xy_int_profile(tile_image, ...
        dt_max_sub, vsl_r_est, ori_vec, voxel_size);
    
    int_profile_str_z = fun_radius_estimation_get_z_int_profile(tile_image, ...
        dt_max_sub, vsl_r_est, ori_vec(3), int_profile_str_xy.min_int, ...
        voxel_size);
    %% Generate intensity profiles for each PSF
    % Find the closest simulated vessel parameters
    est_psf_str = fun_radius_estimation_get_best_fit_psf(psf_est_int, vsl_r_est, ...
        ori_vec_z, int_profile_str_xy, int_profile_str_z);
    assert(isfinite(est_psf_str.corr_max), 'Angle-weighted correlation coefficient is not finite');
    if est_psf_str.corr_max > best_fit_corr
        best_fit_corr = est_psf_str.corr_max;   
        best_fit_r = vsl_r_est;
        best_fit_psf_str = est_psf_str;        
        est_psf_int_interpolation = est_psf_str.best_fit_interpolation;
    else
%         fprintf('New estimation does not increase the fitting accuracy. Return\n');
%         if ori_vec_z < 1/2 && best_fit_corr > 0.99 && iter_est >= 2
%             fprintf('Visualization\n');
%         end
        break;
    end
    if verbose_Q
        fprintf('Best fit PSF FWHM (%.1f, %.1f, %.1f).\n', est_psf_str.best_fit_FWHM);
        fprintf('Lateral intensity profile correlation %.2f\n', est_psf_str.xy_corr_max);
        fprintf('Axial intensity profile correlation %.2f\n', est_psf_str.z_corr_max);
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
        ax_hdl_1.XTick = 1 : vis_tick_step : numel(psf_est_int.psf_FWHM_z_list);
        ax_hdl_1.XTickLabel = arrayfun(@(x) num2str(x, '%.2f'), psf_est_int.psf_FWHM_z_list(1:vis_tick_step:end), ...
            'UniformOutput', false);
        ax_hdl_1.YTick = 1 : vis_tick_step : numel(psf_est_int.psf_FWHM_xy_list);
        ax_hdl_1.YTickLabel = arrayfun(@(x) num2str(x, '%.2f'), psf_est_int.psf_FWHM_xy_list(1:vis_tick_step:end), ...
            'UniformOutput', false);
        ax_hdl_1.XLabel.String = 'FWHM_z(\mum)';
        ax_hdl_1.YLabel.String = 'FWHM_{xy}(\mum)';
        ax_hdl_1.Title.String = sprintf('Avg. Corr. (max = %.5f)', best_fit_psf_str.corr_max);
        
        ax_hdl_2 = subplot(2,3,2);
        imagesc(ax_hdl_2, best_fit_psf_str.xy_corr_mat);
        colormap(ax_hdl_2, 'jet');
        cbar_hdl = colorbar(ax_hdl_2);
        cbar_hdl.Label.String = 'Correlation';
        ax_hdl_2.XTick = 1 : vis_tick_step : numel(psf_est_int.psf_FWHM_z_list);
        ax_hdl_2.XTickLabel = arrayfun(@(x) num2str(x, '%.2f'), psf_est_int.psf_FWHM_z_list(1:vis_tick_step:end), ...
            'UniformOutput', false);
        ax_hdl_2.YTick = 1 : vis_tick_step : numel(psf_est_int.psf_FWHM_xy_list);
        ax_hdl_2.YTickLabel = arrayfun(@(x) num2str(x, '%.2f'), psf_est_int.psf_FWHM_xy_list(1:vis_tick_step:end), ...
            'UniformOutput', false);
        ax_hdl_2.XLabel.String = 'FWHM_z(\mum)';
        ax_hdl_2.YLabel.String = 'FWHM_{xy}(\mum)';
        ax_hdl_2.Title.String = 'Corr. XY';
        
        ax_hdl_3 = subplot(2,3,3);
        imagesc(ax_hdl_3, best_fit_psf_str.z_corr_mat);
        colormap(ax_hdl_3, 'jet');
        cbar_hdl = colorbar(ax_hdl_3);
        cbar_hdl.Label.String = 'Correlation';
        ax_hdl_3.XTick = 1 : vis_tick_step : numel(psf_est_int.psf_FWHM_z_list);
        ax_hdl_3.XTickLabel = arrayfun(@(x) num2str(x, '%.2f'), psf_est_int.psf_FWHM_z_list(1:vis_tick_step:end), ...
            'UniformOutput', false);
        ax_hdl_3.YTick = 1 : vis_tick_step : numel(psf_est_int.psf_FWHM_xy_list);
        ax_hdl_3.YTickLabel = arrayfun(@(x) num2str(x, '%.2f'), psf_est_int.psf_FWHM_xy_list(1:vis_tick_step:end), ...
            'UniformOutput', false);
        ax_hdl_3.XLabel.String = 'FWHM_z(\mum)';
        ax_hdl_3.YLabel.String = 'FWHM_{xy}(\mum)';
        ax_hdl_3.Title.String = sprintf('Corr. Z (v_z = %.2f)', ori_vec_z);
        
        ax_hdl_4 = subplot(2,3,4);
        plot(ax_hdl_4, int_profile_str_xy.r, int_profile_str_xy.int_val_n_smooth);
        hold(ax_hdl_4, 'on');
        plot(ax_hdl_4, best_fit_psf_str.xy_max_corr_int_r + best_fit_psf_str.xy_max_corr_shift_um, best_fit_psf_str.xy_max_corr_int_profile);
        ax_hdl_4.XLabel.String = 'r (\mum)';
        ax_hdl_4.YLabel.String = 'Normalized intensity';
        legend('Measurement', 'Best fit');
        ax_hdl_4.Title.String = sprintf('Lateral intensity profile\nr_{est} = %.1f \\mum, PSF FWHM_{xy} = %.1f \\mum', ...
            vsl_r_est, best_fit_psf_str.best_fit_FWHM(1));
        
        ax_hdl_5 = subplot(2,3,5);
        plot(ax_hdl_5, int_profile_str_z.r, int_profile_str_z.int_val_n_smooth);
        hold(ax_hdl_5, 'on');
        plot(ax_hdl_5, best_fit_psf_str.z_max_corr_int_r + best_fit_psf_str.z_max_corr_shift_um, best_fit_psf_str.z_max_corr_int_profile);
        ax_hdl_5.XLabel.String = 'r (\mum)';
        ax_hdl_5.YLabel.String = 'Normalized intensity';
        legend('Measurement', 'Best fit');
        ax_hdl_5.Title.String = sprintf('Axial intensity profile\nr_{est} = %.1f \\mum, PSF FWHM_z = %.1f \\mum', ...
            vsl_r_est, best_fit_psf_str.best_fit_FWHM(3));     
        %% Save visualization
        dataset_name = 'WholeBrain';
        stack = 'ML20200201';
        DataManager = FileManager;
        im_save_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), ...
            'PSF_estimation');
        fig_fp = fullfile(im_save_folder, sprintf('Corr_map_and_best_fit_int_profile_vsl_r_est_%d_nm_z_comp_%d.png', ...
            round(1000 * vsl_r_est), round(ori_vec_z * 100)));
        fun_print_image_in_several_formats(fig_hdl, fig_fp);
    end
end
est_psf_str = best_fit_psf_str;
est_psf_str.vsl_r_est = best_fit_r;
est_psf_str.vsl_ori_vec = ori_vec;
end