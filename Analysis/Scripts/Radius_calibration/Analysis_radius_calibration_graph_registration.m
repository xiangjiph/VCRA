clc;clear;close all;
DataManager = FileManager;
%% Parameters
registration_method = 'rigid';
skel_reg_opt.method = registration_method;
skel_reg_opt.rot = 1;
skel_reg_opt.scale = 1;
skel_reg_opt.viz = false;              % show every iteration
skel_reg_opt.outliers = 0.7;       % use 0.7 noise weight
skel_reg_opt.max_it = 100;        % Maxinum number of iteration
skel_reg_opt.fgt = 0;              % do not use FGT (default)
skel_reg_opt.tol = 1e-8;         % Error tolorance ( how is it defined? )
skel_reg_opt.normalize = 1;        % normalize to unit variance and zero mean before registering (default)
skel_reg_opt.corresp = 1;          % compute correspondence vector at the end of registration (not being estimated by default)
skel_reg_opt.max_num_desc = 2e4;
skel_reg_opt.rm_descriptor_by_pdist2_Q = false;
skel_reg_opt.pixshift_search_shift_z = [20, 40, 60];
skel_reg_opt.scan_pixshift_Q = true;
skel_reg_opt.only_points_in_overlap_vol_Q = false;

skel_reg_opt.establish_correspondance_method = 'minimum_distance';
skel_reg_opt.brute_search_min_med_dist = 5;
skel_reg_opt.min_cc_inconsistancy = 0.1;
skel_reg_opt.max_pair_dist = 10;
% skel_reg_opt.establish_correspondance_method = 'mutual_closest';
min_depth_um = 25;
max_depth_um = 400;
skel_pre_downsampling_step = 2;
% FFT registration
vis_FFT_Q = false;
mask_fft_ds_ratio = 0.5;
% Parallel computing:
num_processes = 16;
p_hdl = gcp('nocreate');
delete(p_hdl);
p_hdl = parpool(num_processes);
num_thread_per_prop = ceil(48  / num_processes);
%%
dataset_name = 'Vessel_radius_calibration';
stack_list = {'DK20200504_WT2', 'DK20200504_WT1', 'DK20200511_WT5', 'DK20200517_WT7', 'DK20200427_WT3', 'DK20200427_WT4'};
% stack_list = {'DK20200504_WT2', 'DK20200511_WT5', 'DK20200517_WT7', 'DK20200504_WT1'};
num_stack = numel(stack_list);

for iter_dataset = 1 : num_stack
    stack = stack_list{iter_dataset};
    image_group_list = {'In_vivo', 'Post_perfusion'};
    %% Match ROI - assume roi index does not match
    image_group_info_cell = cell(size(image_group_list));
    for iter_group = 1 : numel(image_group_list)
        image_group = image_group_list{iter_group};
        image_group_info_cell{iter_group} = DataManager.load_data(fullfile(...
            DataManager.fp_raw_data_folder(dataset_name, stack), image_group, ...
            sprintf('image_group_info.mat')));
    end
    assert(image_group_info_cell{1}.num_stack == image_group_info_cell{2}.num_stack);
    num_roi = image_group_info_cell{1}.num_stack ;
    registration_result = cell(num_roi, 1);
    mov_group = image_group_list{2};
    fix_group = image_group_list{1};
    
    mov_grid_template = sprintf('240_cube_%s_%%d', mov_group);
    fix_grid_template = sprintf('240_cube_%s_%%d_re', fix_group);
    %% Parallel - skeleton registration
    parfor iter_roi = 1 : num_roi
        maxNumCompThreads(num_thread_per_prop);
        try
            grid_version_fixed = sprintf(fix_grid_template, iter_roi);
            grid_version_mov = sprintf(mov_grid_template, iter_roi);
            vg_fix = DataManager.load_graph_in_block(dataset_name, stack, grid_version_fixed, 0, 0, 0);
            vg_mov = DataManager.load_graph_in_block(dataset_name, stack, grid_version_mov, 0, 0, 0);
            %% Initialize displacement vector
            if isfield(vg_fix.info, 'data_info')
                disp_vec_0 = [vg_fix.info.data_info.slice1_11_disp_vec_wrt_im_11_yx; 0].' .* vg_fix.info.voxel_size_um;
            elseif isfield(vg_mov.info, 'data_info')
                disp_vec_0 = - [vg_mov.info.data_info.slice1_11_disp_vec_wrt_im_11_yx; 0].' .* vg_mov.info.voxel_size_um;
            end
            %% Direct point cloud registration
            fix_sub = fun_ind2sub(vg_fix.num.mask_size, vg_fix.link.pos_ind) .* vg_fix.info.voxel_size_um;
            mov_sub = fun_ind2sub(vg_mov.num.mask_size, vg_mov.link.pos_ind) .* vg_mov.info.voxel_size_um;
            
            fix_selected_idx = find(fix_sub(:, 3) >= min_depth_um & ...
                fix_sub(:, 3) <= max_depth_um);
            mov_selected_idx = find(mov_sub(:, 3) >= min_depth_um & ...
                mov_sub(:, 3) <= max_depth_um);
            %             invivo_selected_idx = 1 : size(invivo_sub, 1);
            %             perfusion_selected_idx = 1 : size(perfusion_sub, 1);
            registered_str = fun_match_vessel_skeleton(...
                fix_sub(fix_selected_idx, :), mov_sub(mov_selected_idx, :), ...
                vg_fix.link.label(fix_selected_idx), vg_mov.link.label(mov_selected_idx), ...
                disp_vec_0, 3, skel_reg_opt);
            
            % To do: Use registration result to initialize large vessel
            % resgistration       
            registered_str.Fixed_matched_input_idx = fix_selected_idx(registered_str.Fixed_matched_input_idx);
            registered_str.Moving_matched_input_idx = mov_selected_idx(registered_str.Moving_matched_input_idx);
            
            assert(all(registered_str.Fixed_cc_label == vg_fix.link.label(registered_str.Fixed_matched_input_idx)), ...
                'Output cc labels are inconsistent with the input labels');
            assert(all(registered_str.Moving_cc_label == vg_mov.link.label(registered_str.Moving_matched_input_idx)), ...
                'Output cc labels are inconsistent with the input labels');
            %%
            registered_str = fun_analysis_post_process_registered_vessel_graph(...
                registered_str, vg_fix, vg_mov);
            DataManager.write_data(registered_str.filepath, registered_str);
            % Save for visualization
            registration_result{iter_roi} = registered_str;
        catch ME
            fprintf('Fail to align ROI %d\n', iter_roi);
            fprintf('Error message: %s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
            fprintf('Skip the following script.\n');
            continue;
        end
    end
    
    %% Visualization
    for iter_roi = 1 : num_roi
        registered_str = registration_result{iter_roi};
        if ~isempty(registered_str.Fixed_sub)
            fig_hdl = figure;
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) *  2.5;
            tiledlayout(fig_hdl, 2, 2);
            ax_hdl = nexttile;
            scatter3(ax_hdl, registered_str.Fixed_sub(:, 1), registered_str.Fixed_sub(:, 2), registered_str.Fixed_sub(:, 3), '.');
            hold(ax_hdl, 'on');
            scatter3(ax_hdl, registered_str.Moving_sub(:, 1), registered_str.Moving_sub(:, 2), registered_str.Moving_sub(:, 3), 'o');
            ax_hdl.DataAspectRatio = [1,1,1];
            ax_hdl.XLim(1) = 0;
            ax_hdl.YLim(1) = 0;
            ax_hdl.ZLim(1) = 0;
            ax_hdl.ZDir = 'reverse';
            ax_hdl.XLabel.String = 'X (\mum)';
            ax_hdl.YLabel.String = 'Y (\mum)';
            ax_hdl.ZLabel.String = 'Z (\mum)';
            ax_hdl.FontSize = 14;
            leg_hdl = legend(ax_hdl, strrep(fix_group, '_', ' '), strrep(mov_group, '_', ' ' ));
            leg_hdl.Location = 'northeast';
            
            ax_hdl_2 = nexttile;
            histogram(ax_hdl_2, registered_str.Dist_fixed_2_moving, 0 : 0.5 : 10, 'Normalization', 'count');
            leg_string = sprintf('Number of pairs: %d\nRigid scale factor %.4f\nMean displacement: %.2f (\\mum)\nMedian displacement %.2f (\\mum)', ...
                size(registered_str.Dist_fixed_2_moving, 1), registered_str.s, mean(registered_str.Dist_fixed_2_moving), median(registered_str.Dist_fixed_2_moving));
            leg_hdl_2 = legend(ax_hdl_2, leg_string, 'Location', 'northeast');
            ax_hdl_2.XLabel.String = 'Residual displacement (\mum)';
            ax_hdl_2.YLabel.String = 'Counts';
            ax_hdl_2.FontSize = 14;
            box(ax_hdl_2, 'off');
            
            ax_hdl_3 = nexttile;
            vis_r_idx = 3;
            hist_r_edge = [0 : 0.2 : 3, 3.5 : 0.5 : 5];
            histogram2(ax_hdl_3, registered_str.Fixed_r_um(:, vis_r_idx), registered_str.Moving_r_um, hist_r_edge, hist_r_edge, 'DisplayStyle', 'tile');
            cbar_hdl = colorbar(ax_hdl_3);
            cbar_hdl.Label.String = 'Number of data points';
            ax_hdl_3.XLabel.String = sprintf('%s radius (\\mum)', strrep(fix_group, '_', ' '));
            ax_hdl_3.YLabel.String = sprintf('%s radius (\\mum)', strrep(mov_group, '_', ' '));
            ax_hdl_3.FontSize = 14;
            ax_hdl_3.Title.String = 'Matched voxel radius comparison';
            ax_hdl_3.DataAspectRatio = [1,1,1];
            % Linear regression:
            % fit_hdl = fitlm(X_r, Y_r, 'Intercept', false);
            
            % Matched vessel segment
            ax_hdl_4 = nexttile;
            histogram2(ax_hdl_4, registered_str.Fixed_cc_features.radius_median(:, vis_r_idx), ...
                registered_str.Moving_cc_features.radius_median, hist_r_edge, hist_r_edge, 'DisplayStyle', 'tile');
            cbar_hdl = colorbar(ax_hdl_4);
            cbar_hdl.Label.String = 'Number of data points';
            ax_hdl_4.XLabel.String = sprintf('%s radius (\\mum)', strrep(fix_group, '_', ' '));
            ax_hdl_4.YLabel.String = sprintf('%s radius (\\mum)', strrep(mov_group, '_', ' '));
            ax_hdl_4.FontSize = 14;
            ax_hdl_4.Title.String = 'Matched segment radius comparison';
            ax_hdl_4.DataAspectRatio = [1,1,1];
            %%
            fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack),...
                'Invivo_vs_postPerfusion_v2', sprintf('%s_%s_%d_regid_%s_to_%s_stat.png',...
                dataset_name, stack, iter_roi, mov_group, fix_group));
            fun_print_image_in_several_formats(fig_hdl, fig_fp);
            delete(fig_hdl);
        end
    end
    fprintf('Finish aligning vessel skeleton for all the images in %s\n', stack);
end