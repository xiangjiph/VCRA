set_env
DataManager = FileManager;
%% Parameters
opt_mask2graph = struct;
% Parameter for skeleton recentering
opt_mask2graph.max_rc_int = 50000;
% Parameter for filling holes in the segmentation mask
opt_mask2graph.node_bbox_expand = 20;

opt_auto_refine = struct;
opt_auto_refine.max_self_loop_length = 30;
opt_auto_refine.internal_offset = 16;
opt_auto_refine.pruning_max_length = 10;
opt_auto_refine.max_bilink_loop_length = 15;

target_isotropic_voxel_size_um = 1;

medfilt_tile_image_Q = true;
%% Load PSF simulation data
% psf_est_int = DataManager.load_data(DataManager.fp_metadata_file('WholeBrain', 'ML_2018_08_15', 'psf_fitting_data'));
% % psf_fp = fullfile(DataManager.fp_metadata_file('DKLab', 'Rui_PSF_simulation', 'psf_fitting_data'));
% % psf_est_int = load(psf_fp);
% psf_FWHM_array = cat(1, psf_est_int.psf_FWHM{:});
% 
% num_psf = numel(psf_est_int.psf_FWHM);
% for iter_psf = 1 : num_psf
%     tmp_psf = psf_est_int.edge_int_str{iter_psf};
%     tmp_num_dist = numel(tmp_psf.normalized_radial_int_dist);
%     tmp_psf.normalized_radial_int_dist_itp = cell(size(tmp_psf.normalized_radial_int_dist));
%     tmp_psf.normalized_vertical_radial_int_dist_itp = cell(size(tmp_psf.normalized_vertical_radial_int_dist));
%     for iter_dist = 1 : tmp_num_dist
%         tmp_int_xy_dist = tmp_psf.normalized_radial_int_dist{iter_dist};
%         tmp_psf.normalized_radial_int_dist_itp{iter_dist} = griddedInterpolant(tmp_int_xy_dist(:, 1), ...
%             tmp_int_xy_dist(:, 2), 'linear', 'nearest');
%         
%         tmp_int_z_dist = tmp_psf.normalized_vertical_radial_int_dist{iter_dist};
%         tmp_psf.normalized_vertical_radial_int_dist_itp{iter_dist} = griddedInterpolant(tmp_int_z_dist(:, 1), ...
%             tmp_int_z_dist(:, 2), 'linear', 'nearest');
%     end
%     psf_est_int.edge_int_str{iter_psf} = tmp_psf;
% end
% fprintf('Finish pre-computing the intensity interpolation\n');
%%
initial_est_PSF_FWHM_xy = 0.9;
initial_est_PSF_FWHM_z = 5.7;
simu_dataset_name = 'Vessel_radius_calibration';
simu_stack = 'Post_perfusion_simulation';
simu_int_fp = DataManager.fp_metadata_file(simu_dataset_name, ...
    simu_stack, 'Estimated_PSF_edge_intensity');
if ~isfile(simu_int_fp)
    gpuDevice(2);
    vessel_radius_list = [0.3 : 0.1 : 15];
    z_component_list = 0 : 0.05 : 1;
    theta_list = asin(z_component_list);
    
    psf_sigma = [ones(1,2) .* initial_est_PSF_FWHM_xy, initial_est_PSF_FWHM_z] ./ (2 * sqrt(2 * log(2)));
    
    dx_min_scale_ratio = 0.1;
    tmp_tic = tic;
    est_psf_str = fun_radius_estimation_compute_vessel_edge_int_given_psf(psf_sigma, ...
        vessel_radius_list, theta_list, dx_min_scale_ratio, false);
    est_psf_str.vsl_ori_vec_z_list = z_component_list;
    est_psf_str.dataset_name = simu_dataset_name;
    est_psf_str.stack = simu_stack;
    est_psf_str.filepath = simu_int_fp;
    DataManager.write_data(est_psf_str.filepath, est_psf_str);
    fprintf('Finish writing the simulated intensity profiles to %s\n', ...
        est_psf_str.filepath);
    toc(tmp_tic);
else
    est_psf_str = DataManager.load_data(simu_int_fp);
end
est_psf_str = est_psf_str.n_min_edge_int_interpolation;
% fprintf('Finish computation. Elapsed time is %f seconds.\n', toc(tmp_tic));

%% Local two photon image peprocessing
% Convert raw iamge to frame-averaged image. Correct pixel shift due to scanner mechanical recovery error.
dataset_name = 'Vessel_radius_calibration';
average_method = 'Average';
image_group_list = {'Post_perfusion'};
stack_list = {'DK20200504_WT2', 'DK20200427_WT3', 'DK20200427_WT4', 'DK20200504_WT1', 'DK20200511_WT5', 'DK20200517_WT7'};
num_image_group = numel(image_group_list);
num_data_stack = numel(stack_list);
%%
num_processes = 6;
p_hdl = gcp('nocreate');
delete(p_hdl);
p_hdl = parpool(num_processes);
num_thread_per_prop = 8;

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
        parfor iter_roi = 1 : num_roi
            use_gpu_id = double(mod(iter_roi, num_processes) > 2) + 1;
            maxNumCompThreads(num_thread_per_prop);
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
            %% Image enhancement
            tmp_tic = tic;
            avg_data_enhance = fun_stretch_contrast(medfilt3(avg_data), 0.005, 0.999);
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
            tmp_grid_version = sprintf('240_cube_%s_%d', image_group, iter_roi);
            tmp_grid_info.version = tmp_grid_version;
            tmp_grid_info.data_type = im_data_type;
            tmp_grid_info.voxel_size_um = tmp_rz_voxel_size_um;
            DataManager.write_grid_info(tmp_grid_info, tmp_grid_info.dataset_name, ...
                tmp_grid_info.stack, tmp_grid_info.version);
            
            exit_code = fun_generate_block_data_from_image_stack(tmp_avg_im_rz, tmp_grid_info);
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
            task_str.gpuDevice = use_gpu_id;
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
            seg_parameters.morp_min_cc_size = 64;
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
            exit_code = fun_task_segmentation(task_str);
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
            vessel_graph.info.grid_c_vergion = tmp_grid_version;
            vessel_graph.info.grid_version = tmp_grid_version;
            vessel_graph.info.voxel_size_um = tmp_rz_voxel_size_um;
            vessel_graph.info.raw_data_fp = tmp_filename;
            vessel_graph.info.image_group = image_group;
            vessel_graph.info.ROI_ID = iter_roi;
            vessel_graph.info.RE_mefilt3_Q = medfilt_tile_image_Q;
            vessel_graph.info.segmentation_parameter = seg_parameters;
            %% PSF - Radius joint estimation
            % Due to the computational cost, only do for part of the vessels?
            % The image shouldn't be saturated. So, use the original averaged image
%             if medfilt_tile_image_Q
%                 tile_image = medfilt3(avg_data);
%             else
%                 tile_image = avg_data;
%             end
%             link_cc_sub = cellfun(@(x) fun_ind2sub(vessel_graph.num.mask_size, ...
%                 x), vessel_graph.link.cc_ind, 'UniformOutput', false);
%             link_cc_sub_COM = cellfun(@(x) mean(x, 1), link_cc_sub, 'UniformOutput', false);
%             link_cc_sub_COM = cat(1, link_cc_sub_COM{:});
%             selected_link_for_PSF_est_label = find(link_cc_sub_COM(:, 3) < 800);
%             %     assert(numel(selected_link_for_PSF_est_label) == numel(unique(selected_link_for_PSF_est_label)));
%             skeleton_cc_sub_um = cellfun(@(x) fun_ind2sub(vessel_graph.num.mask_size, x) .* ...
%                 vessel_graph.info.voxel_size_um, vessel_graph.link.cc_ind(selected_link_for_PSF_est_label), 'UniformOutput', false);
%             skeleton_cc_r_um = cellfun(@(x) double(full(vessel_graph.radius(x))), ...
%                 vessel_graph.link.cc_ind(selected_link_for_PSF_est_label), 'UniformOutput', false);
%             psf_stat = fun_radius_estimation_get_PSF_stat_in_image_stack(psf_est_int, ...
%                 tile_image, voxel_size_um, skeleton_cc_sub_um, skeleton_cc_r_um);
            %         %
            %         %     psf_stat_medfilt = fun_radius_estimation_get_PSF_stat_in_image_stack(psf_est_int, ...
            %         %         medfilt3(tile_image), voxel_size_um, skeleton_cc_sub_um, skeleton_cc_r_um);
                    %% Visualize PSF estimation result
%                         tmp_vis_data = psf_stat;
%                         tmp_vis_selected_Q = link_cc_sub_COM(:, 3) < 250;
%                         vis_link_ori_vec = cat(1, tmp_vis_data.link_ori_vec{tmp_vis_selected_Q});
%                         vis_best_fit_FWHM = cat(1, tmp_vis_data.best_fit_PSF_FWHM{tmp_vis_selected_Q});
%                         vis_r = cat(1, tmp_vis_data.refine_r_est{tmp_vis_selected_Q});
%                         vis_best_fit_corr = cat(1, tmp_vis_data.best_fit_corr{tmp_vis_selected_Q});
%             
%                         tmp_vis_Q = (vis_r <= 4 & vis_best_fit_corr > 0.99 & ...
%                             abs(vis_link_ori_vec(:, 3)) < sqrt(1/2));
%                         % tmp_vis_Q = ~isnan(best_fit_corr);
%                         fig_hdl = figure;
%                         fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 1] .* 1.5;
%                         ax_hdl_1 = subplot(1,3,1);
%                         tmp_data = vis_best_fit_FWHM(tmp_vis_Q, 1);
%                         histogram(ax_hdl_1, tmp_data, psf_est_int.psf_FWHM_xy_list - 0.05, 'Normalization', 'pdf');
%                         legend(ax_hdl_1, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data)));
%                         ax_hdl_1.XLabel.String = 'FWHM_{xy} (\mum)';
%                         ax_hdl_1.YLabel.String = 'PDF';
%                         ax_hdl_1.FontSize = 14;
%                         ax_hdl_2 = subplot(1,3,2);
%                         tmp_data = vis_best_fit_FWHM(tmp_vis_Q, 3);
%                         histogram(ax_hdl_2, tmp_data, [0, psf_est_int.psf_FWHM_z_list], 'Normalization', 'pdf');
%                         legend(ax_hdl_2, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data)), 'Location', 'northwest');
%                         ax_hdl_2.XLabel.String = 'FWHM_{z} (\mum)';
%                         ax_hdl_2.YLabel.String = 'PDF';
%                         ax_hdl_2.FontSize = 14;
%             
%                         ax_hdl_3 = subplot(1,3,3);
%                         histogram(ax_hdl_3, vis_r(tmp_vis_Q), 'Normalization', 'pdf');
%                         legend(ax_hdl_3, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(vis_r(tmp_vis_Q))), 'Location', 'northeast');
%                         ax_hdl_3.XLabel.String = 'Estimated radius (\mum)';
%                         ax_hdl_3.YLabel.String = 'PDF';
%                         ax_hdl_3.FontSize = 14;
%             %
%                         if medfilt_tile_image_Q
%                             fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'PSF_estimation', image_group, sprintf('%s_%s_%s_%s_%d_PSF_FWHM_estimation_hists.png', ...
%                                 dataset_name, stack, image_group, average_method, iter_roi));
%                         else
%                             fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'PSF_estimation_medfilt', image_group, sprintf('%s_%s_%s_%s_%d_PSF_FWHM_estimation_hists.png', ...
%                                 dataset_name, stack, image_group, average_method, iter_roi));
%                         end
%                         fun_print_image_in_several_formats(fig_hdl, fig_fp);
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
            est_radius_str = fun_radius_estimation_in_image_stack(est_psf_str, ...
                tile_image, voxel_size_um, skeleton_cc_sub_um, skeleton_cc_r_um);
            vessel_graph.link.features.radius_median = est_radius_str.cc_stat.median;
            vessel_graph.link.features.radius_mean = est_radius_str.cc_stat.mean;
            vessel_graph.link.features.radius_std = est_radius_str.cc_stat.std;
            
            vessel_graph.link.features.cc_com_sub = cellfun(@(x) mean(x, 1), ...
                skeleton_cc_sub_um, 'UniformOutput', false);
            vessel_graph.link.features.cc_com_sub = cat(1, vessel_graph.link.features.cc_com_sub{:});
            vessel_graph.link.cc_r = est_radius_str.r;
            vessel_graph.link.cc_int = cellfun(@(x) vessel_image(x), vessel_graph.link.cc_ind, 'UniformOutput', false);
            vessel_graph.link.cc_r_valid_Q = est_radius_str.valid_Q;
            vessel_graph.link.cc_local_ori_vec_z = est_radius_str.local_ori_vec_z;
            %% Save vessel graph
            DataManager.write_graph_in_block(vessel_graph, dataset_name, stack, ...
                tmp_grid_version, 0, 0, 0);
            fprintf('Finish converting segmentation to graph. Elapsed time is %f secodns\n', toc(tmp_tic));
            %% Visualize image
            %     fig_hdl = figure;
            %     ax_hdl = axes(fig_hdl);
            %     for iter_sec = 1 : size(tile_image, 3)
            %         iter_sec = 300;
            %         imagesc(ax_hdl, tile_image(:, :, iter_sec));
            %         ax_hdl.DataAspectRatio = [1,1,1];
            %         colorbar(ax_hdl);
            % %         colormap('jet');
            %         pause(0.1);
            %     end
        end
    end
end