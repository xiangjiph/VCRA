set_env
DataManager = FileManager;
%% Parameters
target_isotropic_voxel_size_um = 1;

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

num_processes = 6;
p_hdl = gcp('nocreate');
delete(p_hdl);
p_hdl = parpool(num_processes);

num_thread_per_prop = 8;
%% Local two photon image peprocessing
% Convert raw iamge to frame-averaged image. Correct pixel shift due to scanner mechanical recovery error.
dataset_name = 'Vessel_radius_calibration';
average_method = 'Rank_3_w_gmc';
% average_method = 'Average';
image_group_list = {'In_vivo'};
stack_list = {'DK20200504_WT2', 'DK20200504_WT1', 'DK20200511_WT5', 'DK20200517_WT7', 'DK20200427_WT3', 'DK20200427_WT4'};
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
        parfor iter_roi = 1 : num_roi
            use_gpu_id = double(mod(iter_roi, num_processes) > 2) + 1;
            maxNumCompThreads(num_thread_per_prop);
            try
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
                % Image enhancement is necessary. The Tiff images were
                % stored in uint16, but the acquisition was actually int16.
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
                task_str.gpuDevice = use_gpu_id;
                [~, task_str] = fun_task_get_task_str(task_str, 1, true);
                
                seg_parameters = struct;
                seg_parameters.voxel_length_um = round(mean(tmp_grid_info.voxel_size_um));
                seg_parameters.rod_filter_radius_um = 1;
                seg_parameters.rod_filter_length_um = round(6*seg_parameters.rod_filter_radius_um + 1);
                seg_parameters.rod_filter_num_omega = 6;
                seg_parameters.vesselness.DoG_scale_list_um = [0.5, 1, 2]./ seg_parameters.voxel_length_um;
                seg_parameters.vesselness_th = 0.1;
                seg_parameters.adp_th_scale_1_um = 8;
                seg_parameters.adp_th_scale_2_um = 16;
                seg_parameters.morp_min_cc_size = 27;
                seg_parameters.max_pool_size = 8;
                seg_parameters.min_bg_std = 250;
                seg_parameters.grid_info = tmp_grid_info;
                seg_parameters.voxel_length_um = round(mean(tmp_grid_info.voxel_size_um));
                seg_parameters.data_type = tmp_grid_info.data_type;
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
                vessel_graph.info.grid_version = tmp_grid_version;
                vessel_graph.info.grid_info = tmp_grid_info;
                vessel_graph.info.voxel_size_um = tmp_rz_voxel_size_um;
                vessel_graph.info.raw_data_fp = tmp_filename;
                vessel_graph.info.image_group = image_group;
                vessel_graph.info.ROI_ID = iter_roi;
                vessel_graph.info.segmentation_parameter = seg_parameters;
                vessel_graph.info.data_info = avg_data_info;
                fprintf('Finish converting segmentation to vessel graph\n');
                %% Save vessel graph
                DataManager.write_graph_in_block(vessel_graph, dataset_name, stack, ...
                    tmp_grid_version, 0, 0, 0);
                fprintf('Finish converting segmentation to graph. Elapsed time is %f secodns\n', toc(tmp_tic));
            catch ME
                fprintf('Unable to process ROI %d.\n', iter_roi);
                rethrow(ME);
            end
        end
    end
    fprintf('Finish processing %s %s\n', dataset_name, stack);
end