set_env
DataManager = FileManager;
%% Parameters
opt_mask2graph = struct;
% Parameter for skeleton recentering
opt_mask2graph.max_rc_int = 0;
% Parameter for filling holes in the segmentation mask
opt_mask2graph.node_bbox_expand = 20;

opt_auto_refine = struct;
opt_auto_refine.max_self_loop_length = 30;
opt_auto_refine.internal_offset = 16;
opt_auto_refine.pruning_max_length = 10;
opt_auto_refine.max_bilink_loop_length = 15;

target_isotropic_voxel_size_um = 1;

medfilt_tile_image_Q = true;
%%
simulation_file_fp = fullfile(DataManager.fp_metadata_file('Vessel_radius_calibration', 'In_vivo_image_simulation', 'In_vivo_image_simulation'));
invivo_im_prof_str = DataManager.load_data(simulation_file_fp);

psf_est_int = DataManager.load_data(DataManager.fp_metadata_file('WholeBrain', 'ML_2018_08_15', 'psf_fitting_data'));
initial_est_PSF_FWHM_xy = 0.9;
initial_est_PSF_FWHM_z = 5.7;
[~, initial_est_PSF_FWHM_xy_ind] = min(abs(psf_est_int.psf_FWHM_xy_list - initial_est_PSF_FWHM_xy));
[~, initial_est_PSF_FWHM_z_ind] = min(abs(psf_est_int.psf_FWHM_z_list - initial_est_PSF_FWHM_z));
initial_est_PSF_list_ind = sub2ind(size(psf_est_int.psf_FWHM), initial_est_PSF_FWHM_xy_ind, initial_est_PSF_FWHM_z_ind);
est_psf_str = psf_est_int.edge_int_str{initial_est_PSF_list_ind};
est_psf_str = est_psf_str.n_min_edge_int_interpolation;
%% Preprocess simulation result
num_parameter_combination = numel(invivo_im_prof_str.edge_int_str);
for iter_p = 1 : num_parameter_combination
    tmp_psf = invivo_im_prof_str.edge_int_str{iter_p};
    tmp_num_dist = numel(tmp_psf.normalized_radial_int_dist);
    tmp_psf.normalized_radial_int_dist_itp = cell(size(tmp_psf.normalized_radial_int_dist));
    tmp_psf.normalized_vertical_radial_int_dist_itp = cell(size(tmp_psf.normalized_vertical_radial_int_dist));
    
    for iter_dist = 1 : tmp_num_dist
        tmp_int_xy_dist = tmp_psf.max_normalized_radial_int_dist{iter_dist};
        if ~isempty(tmp_int_xy_dist)
            tmp_psf.normalized_radial_int_dist_itp{iter_dist} = griddedInterpolant(tmp_int_xy_dist(:, 1), ...
                tmp_int_xy_dist(:, 2), 'linear', 'nearest');
            
            tmp_int_z_dist = tmp_psf.max_normalized_vertical_radial_int_dist{iter_dist};
            tmp_psf.normalized_vertical_radial_int_dist_itp{iter_dist} = griddedInterpolant(tmp_int_z_dist(:, 1), ...
                tmp_int_z_dist(:, 2), 'linear', 'nearest');
        end
    end
    invivo_im_prof_str.edge_int_str{iter_p} = tmp_psf;
end
fprintf('Finish pre-computing the intensity interpolation\n');
%% Profile analysis
% Is 50% for capillary good? Assume completely no fluorescence signal from
% RBC? - Plot 2 um, 0 degree,
vis_r_um = 2.5;
vis_z = 0;
[~, vis_r_idx] = min(abs(invivo_im_prof_str.vsl_rds_list_um - vis_r_um));
[~, vis_z_idx] = min(abs(invivo_im_prof_str.vsl_ori_z_comp - vis_z));
vis_plot_profs = cell(size(invivo_im_prof_str.edge_int_str));
for iter_profs = 1 : numel(vis_plot_profs)
    tmp_str = invivo_im_prof_str.edge_int_str{iter_profs};
    tmp_mat = tmp_str.max_normalized_radial_int_dist{vis_z_idx, vis_r_idx};
    vis_plot_profs{iter_profs} = tmp_mat;
end
%%
vis_bc_lyr_t_list = [0.2 : 0.2 : 0.8];

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
for iter_vis = 1 : numel(vis_bc_lyr_t_list)
    vis_bc_lyr_t_um = vis_bc_lyr_t_list(iter_vis);
    vis_rbc_int_n = 0.2;
    [~, vis_bc_lyr_idx] = min(abs(invivo_im_prof_str.vsl_bc_lyr_thkn_um_list - ...
        vis_bc_lyr_t_um));
    vis_bc_lyr_t_um = invivo_im_prof_str.vsl_bc_lyr_thkn_um_list(vis_bc_lyr_idx);
    [~, vis_rbc_int_idx] = min(abs(invivo_im_prof_str.rbc_int_n_list - ...
        vis_rbc_int_n));
    vis_rbc_int_n = invivo_im_prof_str.rbc_int_n_list(vis_rbc_int_idx);
    tmp_plot_vecs = vis_plot_profs{vis_rbc_int_idx, vis_bc_lyr_idx};
    
    plot(ax_hdl, tmp_plot_vecs(:, 1), tmp_plot_vecs(:, 2), 'LineWidth', 1);
    hold(ax_hdl, 'on');
end
ax_hdl.XLabel.String = 'Horizontal radial direction (\mum)';
ax_hdl.YLabel.String = 'Normalized intensity';
leg_hdl = legend(ax_hdl, arrayfun(@(x) num2str(x, '%.2f'), vis_bc_lyr_t_list, 'UniformOutput', false));
leg_hdl.Title.String = 'BLT (\mum)';
%% Observations
% 1. Normalized edge intensity can be greater than 0.5 if the RBC tube
% intensity is too weak (e.g. 0). In that case, thresholding the local
% image profile by a single threshold would be a problem - create a hole
% 2. So, it seems that the estimation should be based on correlation of the
% entire intensity profile.
%
%% Local two photon image peprocessing
% Convert raw iamge to frame-averaged image. Correct pixel shift due to scanner mechanical recovery error.
dataset_name = 'Vessel_radius_calibration';
% average_method = 'Rank_3_w_gmc';
average_method = 'Average';
image_group_list = {'In_vivo'};
stack_list = {'DK20200504_WT2', 'DK20200504_WT1', 'DK20200427_WT3', 'DK20200427_WT4', 'DK20200511_WT5', 'DK20200517_WT7'};
% stack_list = {'DK20200427_WT3', 'DK20200427_WT4'};
num_image_group = numel(image_group_list);
num_data_stack = numel(stack_list);
%%
for iter_data_stack = 1 : num_data_stack
    % stack = 'DK20200504_WT2';
    % image_group_list = {'In_vivo'};
    stack = stack_list{iter_data_stack};
    for iter_IG = 1 : num_image_group
        image_group = image_group_list{iter_IG};
        % image_group = 'In_vivo';
        % image_group = 'Post_perfusion';
        averaged_image_folder = fullfile(DataManager.fp_processed_data(dataset_name, ...
            stack), sprintf('%s_image', average_method), image_group);
        image_group_info = DataManager.load_data(fullfile(DataManager.fp_raw_data_folder(dataset_name, stack), ...
            image_group, 'image_group_info.mat'));
        averaged_image_file_str = dir(fullfile(averaged_image_folder, '*.tiff'));
        num_roi = numel(averaged_image_file_str);
        assert(num_roi == image_group_info.num_stack);
        %%
        for iter_roi = 1 : num_roi
            %% Load data
            tmp_tic = tic;
            tmp_filename = fullfile(averaged_image_folder,...
                sprintf('%s_%s_%s_%s_%d.tiff', dataset_name, stack, image_group, average_method, iter_roi));
            avg_data = DataManager.load_single_tiff(tmp_filename);
            
            tmp_filename = fullfile(averaged_image_folder,...
                sprintf('%s_%s_%s_%s_%d_info.mat', dataset_name, stack, image_group, average_method, iter_roi));
            avg_data_info = DataManager.load_data(tmp_filename);
            avg_im_size = size(avg_data);
            voxel_size_um = avg_data_info.voxel_size_um;
            im_data_type = class(avg_data);
            fprintf('Finish loading averaged image stack. Elapsed time is %f seconds\n', toc(tmp_tic));
            tmp_grid_version = sprintf('240_cube_%s_%d', image_group, iter_roi);
            %% Image enhancement
            tmp_tic = tic;
            avg_data_enhance = fun_stretch_contrast(medfilt3(avg_data), 0.005, 0.995);
            fprintf('Finish enhancing the image. Elapsed time is %f seconds\n', toc(tmp_tic));
            %% Downsampling
            % Convert to isotropic voxel size
            tmp_target_im_size = round(avg_im_size .* voxel_size_um ./ target_isotropic_voxel_size_um);
            tmp_rz_voxel_size_um = avg_im_size .* voxel_size_um ./tmp_target_im_size;
            tmp_scale_up_factor = tmp_rz_voxel_size_um ./ voxel_size_um;
            tmp_apperant_voxel_size = voxel_size_um ./ tmp_rz_voxel_size_um;
            tmp_avg_im_rz = imresize3(avg_data_enhance, tmp_target_im_size);
            % Blocking for segmentation
            tmp_grid_info = fun_get_grid_from_mask(tmp_avg_im_rz, 240, 16, tmp_target_im_size, 0, false);
            tmp_grid_info.dataset_name = dataset_name;
            tmp_grid_info.stack = stack;
            
            tmp_grid_info.version = tmp_grid_version;
            tmp_grid_info.data_type = im_data_type;
            tmp_grid_info.voxel_size_um = tmp_rz_voxel_size_um;
            DataManager.write_grid_info(tmp_grid_info, tmp_grid_info.dataset_name, ...
                tmp_grid_info.stack, tmp_grid_info.version);
            %% Segmentation task setting
            task_name = sprintf('Segmentation_%s_%d', image_group, iter_roi);
            task_folder = DataManager.fp_task_folder(dataset_name, stack, task_name);
            if ~isfolder(task_folder)
                mkdir(task_folder);
            end
            task_str = struct;
            task_str.DataManager = FileManager;
            task_str.dataset_name = dataset_name;
            task_str.stack = stack;
            task_str.grid_version = tmp_grid_version;
            task_str.task_option = [];
            task_str.task_function_name = 'fun_task_segmentation';
            task_str.task_name = task_name;
            task_str.fun_handle = @fun_segmentation_image_to_mask;
            % Use single process
            task_str.task_list = 1 : tmp_grid_info.num_valid_cube;
            task_str.gpuDevice = 2;
            [~, task_str] = fun_task_get_task_str(task_str, 1, true);
            
            seg_parameters = struct;
            seg_parameters.voxel_length_um = round(mean(tmp_grid_info.voxel_size_um));
            seg_parameters.data_type = tmp_grid_info.data_type;
            seg_parameters.rod_filter_radius_um = 1;
            seg_parameters.rod_filter_length_um = round(6*seg_parameters.rod_filter_radius_um + 1);
            seg_parameters.rod_filter_num_omega = 6;
            seg_parameters.vesselness.DoG_scale_list_um = [0.5, 1, 2]./seg_parameters.voxel_length_um;
            seg_parameters.vesselness_th = 0.1;
            seg_parameters.adp_th_scale_1_um = 8;
            seg_parameters.adp_th_scale_2_um = 16;
            seg_parameters.morp_min_cc_size = 27;
            seg_parameters.max_pool_size = 8;
            seg_parameters.min_bg_std = 250;
            seg_parameters.grid_info = tmp_grid_info;
            
            task_str.opt = seg_parameters;
            task_str.overwrite_Q = true;
            %% Segmentation
            %             iter_label = 14;
            %             cube_sub = tmp_grid_info.bbox_grid_sub_list(iter_label, :);
            %             im_cube = DataManager.load_block_data(dataset_name, stack, tmp_grid_version, cube_sub(1), ...
            %                 cube_sub(2), cube_sub(3));
            %             im_cube = fun_stretch_contrast(im_cube, 0, 1);
            %             im_mask = fun_mouselight_segmentation_1um_cube(im_cube, seg_parameters);
            exit_code = fun_generate_block_data_from_image_stack(tmp_avg_im_rz, tmp_grid_info);
            exit_code = fun_task_in_vivo_TPV_segmentation(task_str);
            %% Annotation - fix the in vivo large vessel segmentation
            %% Convert segmentation to graph without annotation
            fprintf('Convert segmentation to graph\n');
            tmp_tic = tic;
            vessel_mask = DataManager.load_blocks_files('mask', dataset_name, stack, ...
                tmp_grid_version, 1 : tmp_grid_info.grid_size(1), 1 : tmp_grid_info.grid_size(2), ...
                1 : tmp_grid_info.grid_size(3));
            vessel_image = DataManager.load_blocks_files('image', dataset_name, stack, ...
                tmp_grid_version, 1 : tmp_grid_info.grid_size(1), 1 : tmp_grid_info.grid_size(2), ...
                1 : tmp_grid_info.grid_size(3));
            assert(all(size(vessel_mask) == size(vessel_image)), 'Vessel mask size is different from averaged image stack');
            vessel_mask = bwareaopen(vessel_mask, seg_parameters.morp_min_cc_size);
            vessel_mask = imclose(vessel_mask, strel('sphere', 2));
            [vessel_graph, vessel_mask_dt] = fun_graph_mask_to_graph(vessel_mask, vessel_image, opt_mask2graph);
            [vessel_graph, pruning_info_str_1] = fun_graph_delete_hairs_and_short_loops(vessel_graph, vessel_image, opt_auto_refine);
            vessel_graph = fun_graph_add_radius(vessel_graph, vessel_mask_dt);
            % Save graph
            vessel_graph.info.dataset_name = dataset_name;
            vessel_graph.info.stack = stack;
            veseel_graph.info.grid_c_vergion = tmp_grid_version;
            vessel_graph.info.grid_version = tmp_grid_version;
            vessel_graph.info.voxel_size_um = tmp_rz_voxel_size_um;
            vessel_graph.info.raw_data_fp = tmp_filename;
            vessel_graph.info.image_group = image_group;
            vessel_graph.info.ROI_ID = iter_roi;
            vessel_graph.info.RE_mefilt3_Q = medfilt_tile_image_Q;
            vessel_graph.info.segmentation_parameter = seg_parameters;
            fprintf('Finish converting segmentation to vessel graph\n');
            %% Visualize reconstructed mask
            %             vessel_mask_recon = uint8(fun_graph_reconstruction(vessel_graph));
            %             vessel_mask_recon(vessel_graph.link.pos_ind) = 2;
            %             vessel_mask_recon(vessel_graph.node.pos_ind) = 3;
            %             DataManager.visualize_itksnap(vessel_image, vessel_mask_recon);
            %% Radius estimation with fixed PSF parameters
            if medfilt_tile_image_Q
                tile_image = medfilt3(avg_data);
            else
                tile_image = avg_data;
            end
            skeleton_cc_sub_um = cellfun(@(x) fun_ind2sub(vessel_graph.num.mask_size, x) .* ...
                vessel_graph.info.voxel_size_um, vessel_graph.link.cc_ind, 'UniformOutput', false);
            skeleton_cc_r_um = cellfun(@(x) double(full(vessel_graph.radius(x))), ...
                vessel_graph.link.cc_ind, 'UniformOutput', false);
            % Use vessel with RBC simulation
            est_radius_str = fun_in_vivo_radius_estimation_in_image_stack(invivo_im_prof_str, ...
                tile_image, voxel_size_um, skeleton_cc_sub_um, skeleton_cc_r_um, 'DT');
%             r_est_fit = fun_in_vivo_radius_estimation_in_image_stack(invivo_im_prof_str, ...
%                 tile_image, voxel_size_um, skeleton_cc_sub_um, skeleton_cc_r_um, 'Fit');
%             % Use vessel without RBC simulation
%             est_radius_str = fun_radius_estimation_in_image_stack(est_psf_str, ...
%                 tile_image, voxel_size_um, skeleton_cc_sub_um, skeleton_cc_r_um);
            
            vessel_graph.link.features.radius_median = est_radius_str.cc_stat.median;
            vessel_graph.link.features.radius_mean = est_radius_str.cc_stat.mean;
            vessel_graph.link.features.radius_std = est_radius_str.cc_stat.std;
            
            vessel_graph.link.features.cc_com_sub = cellfun(@(x) mean(x, 1), ...
                skeleton_cc_sub_um, 'UniformOutput', false);
            vessel_graph.link.cc_r = est_radius_str.r;
            vessel_graph.link.cc_int = cellfun(@(x) vessel_image(x), vessel_graph.link.cc_ind, 'UniformOutput', false);
            vessel_graph.link.cc_r_valid_Q = est_radius_str.valid_Q;
            vessel_graph.link.cc_local_ori_vec_z = est_radius_str.local_ori_vec_z;
            vessel_graph.link.cc_est_parms = est_radius_str.best_fit_parameters;
            vessel_graph.link.cc_est_corr = est_radius_str.best_fit_corr;
            %% Save vessel graph
            DataManager.write_graph_in_block(vessel_graph, dataset_name, stack, ...
                tmp_grid_version, 0, 0, 0);
            fprintf('Finish converting segmentation to graph. Elapsed time is %f secodns\n', toc(tmp_tic));
        end
    end
    fprintf('Finish processing %s %s\n', dataset_name, stack);
end
%% Comparing Estimation accoridng to simulation with vs without RBC
% tmp_bin = 0 : 0.25 : 8;
% tmp_x = r_est_dt.cc_stat.median;
% tmp_y = est_radius_str.cc_stat.median;
% tmp_selected_Q = tmp_x < tmp_bin(end) & tmp_y < tmp_bin(end);
% tmp_x = tmp_x(tmp_selected_Q);
% tmp_y = tmp_y(tmp_selected_Q);
% fig_hdl = figure;
% ax_hdl = axes(fig_hdl);
% histogram2(ax_hdl, tmp_x, tmp_y, tmp_bin, tmp_bin, 'DisplayStyle', 'tile');
% ax_hdl.XLabel.String = 'Profile w RBC (\mum)';
% ax_hdl.YLabel.String = 'Profile w/o RBC (\mum)';
% hold(ax_hdl, 'on');
% plot(ax_hdl, tmp_bin, tmp_bin, 'k-.', 'LineWidth', 1);
% ax_hdl.FontSize = 14;
% ax_hdl.DataAspectRatio = [1,1,1];
% cbar_hdl = colorbar(ax_hdl);
% cbar_hdl.Label.String = 'Number of vessels';
% ax_hdl.ColorScale = 'log';
% 
% tmp_fit = fitlm(tmp_x, tmp_y, 'Intercept', false);
% leg_string = sprintf('Slope: %.2f\nR^2-Adj: %.2f', tmp_fit.Coefficients.Estimate(1), ...
%     tmp_fit.Rsquared.Adjusted);
% leg_hdl = legend(ax_hdl, leg_string, 'Location', 'northwest');
% fun_print_image_in_several_formats(fig_hdl, fullfile(DataManager.fp_visualization_folder(...
%     dataset_name, stack), 'Radius_estimation_comparison', sprintf('%s_%s_%s_%d_re_dt_wRBC_vs_woRBC_medfilt3.png', ...
%     dataset_name, stack, image_group, iter_roi)));
%% Debug
% check_cc_idx = find(abs(r_est_dt.cc_stat.median - est_radius_str.cc_stat.median) > 2);
%% 
% tmp_bin = 0 : 0.25 : 6;
% tmp_x = r_est_dt.cc_stat.median;
% tmp_y = r_est_fit.cc_stat.median;
% tmp_num_vxl_per_cc = cellfun(@numel, r_est_dt.r);
% tmp_selected_Q = tmp_x < tmp_bin(end) & tmp_y < tmp_bin(end) & tmp_num_vxl_per_cc > 8;
% tmp_x = tmp_x(tmp_selected_Q);
% tmp_y = tmp_y(tmp_selected_Q);
% 
% fig_hdl = figure;
% ax_hdl = axes(fig_hdl);
% histogram2(ax_hdl, tmp_x, tmp_y, tmp_bin, tmp_bin, 'DisplayStyle', 'tile');
% ax_hdl.XLabel.String = 'Radius estimated by DT (\mum)';
% ax_hdl.YLabel.String = 'Radius estiamted by profile (\mum)';
% hold(ax_hdl, 'on');
% plot(ax_hdl, tmp_bin, tmp_bin, 'k-.', 'LineWidth', 1);
% ax_hdl.FontSize = 14;
% ax_hdl.DataAspectRatio = [1,1,1];
% cbar_hdl = colorbar(ax_hdl);
% cbar_hdl.Label.String = 'Number of vessels';
% tmp_fit = fitlm(tmp_x, tmp_y, 'Intercept', false);
% ax_hdl.ColorScale = 'log';
% leg_string = sprintf('Slope: %.2f\nR^2-Adj: %.2f', tmp_fit.Coefficients.Estimate(1), ...
%     tmp_fit.Rsquared.Adjusted);
% leg_hdl = legend(ax_hdl, leg_string, 'Location', 'northwest');
% fun_print_image_in_several_formats(fig_hdl, fullfile(DataManager.fp_visualization_folder(...
%     dataset_name, stack), 'Radius_estimation_comparison', sprintf('%s_%s_%s_%d_re_dt_vs_profile_medfilt3.png', ...
%     dataset_name, stack, image_group, iter_roi))); 
% est_diff = r_est_dt.cc_stat.median - r_est_fit.cc_stat.median;
% check_cc_idx = find(abs(est_diff) > 2);
%% Visualize PSF estimation result
%             tmp_vis_data = psf_stat;
%             vis_link_ori_vec = cat(1, tmp_vis_data.link_ori_vec{:});
%             vis_best_fit_FWHM = cat(1, tmp_vis_data.best_fit_PSF_FWHM{:});
%             vis_r = cat(1, tmp_vis_data.refine_r_est{:});
%             vis_best_fit_corr = cat(1, tmp_vis_data.best_fit_corr{:});
%
%             tmp_vis_Q = (vis_r <= 4 & vis_best_fit_corr > 0.99 & ...
%                 abs(vis_link_ori_vec(:, 3)) < sqrt(1/2));
%             % tmp_vis_Q = ~isnan(best_fit_corr);
%             fig_hdl = figure;
%             fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 1] .* 1.5;
%             ax_hdl_1 = subplot(1,3,1);
%             tmp_data = vis_best_fit_FWHM(tmp_vis_Q, 1);
%             histogram(ax_hdl_1, tmp_data, psf_est_int.psf_FWHM_xy_list - 0.05, 'Normalization', 'pdf');
%             legend(ax_hdl_1, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data)));
%             ax_hdl_1.XLabel.String = 'FWHM_{xy} (\mum)';
%             ax_hdl_1.YLabel.String = 'PDF';
%             ax_hdl_1.FontSize = 14;
%             ax_hdl_2 = subplot(1,3,2);
%             tmp_data = vis_best_fit_FWHM(tmp_vis_Q, 3);
%             histogram(ax_hdl_2, tmp_data, [0, psf_est_int.psf_FWHM_z_list], 'Normalization', 'pdf');
%             legend(ax_hdl_2, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data)), 'Location', 'northwest');
%             ax_hdl_2.XLabel.String = 'FWHM_{z} (\mum)';
%             ax_hdl_2.YLabel.String = 'PDF';
%             ax_hdl_2.FontSize = 14;
%
%             ax_hdl_3 = subplot(1,3,3);
%             histogram(ax_hdl_3, vis_r(tmp_vis_Q), 'Normalization', 'pdf');
%             legend(ax_hdl_3, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(vis_r(tmp_vis_Q))), 'Location', 'northeast');
%             ax_hdl_3.XLabel.String = 'Estimated radius (\mum)';
%             ax_hdl_3.YLabel.String = 'PDF';
%             ax_hdl_3.FontSize = 14;
%
%             if medfilt_tile_image_Q
%                 fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'PSF_estimation', image_group, sprintf('%s_%s_%s_ROI_%d_PSF_FWHM_estimation_hists.png', ...
%                     dataset_name, stack, image_group, iter_roi));
%             else
%                 fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'PSF_estimation_medfilt', image_group, sprintf('%s_%s_%s_ROI_%d_PSF_FWHM_estimation_hists.png', ...
%                     dataset_name, stack, image_group, iter_roi));
%             end
%             fun_print_image_in_several_formats(fig_hdl, fig_fp);