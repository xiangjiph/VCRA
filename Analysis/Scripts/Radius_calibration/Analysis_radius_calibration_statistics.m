clc;clear;close all;
DataManager = FileManager;
dataset_name = 'Vessel_radius_calibration';
stack_list = {'DK20200504_WT2', 'DK20200511_WT5', 'DK20200517_WT7', 'DK20200504_WT1',...
    'DK20200427_WT3', 'DK20200427_WT4'};
% stack_list = {'DK20200504_WT2', 'DK20200504_WT1', 'DK20200511_WT5', 'DK20200517_WT7'};
num_stack = numel(stack_list);
stack_data = cell(num_stack, 1);
%%
for iter_stack = 1 : num_stack
    stack = stack_list{iter_stack};
    image_group_list = {'In_vivo', 'Post_perfusion'};
    im_group_mov = image_group_list{1};
    im_group_fixed = image_group_list{2};
    registration_method = 'rigid';
    visualization_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), ...
        'stack_statistics');
    num_image_group = numel(image_group_list);
    
    mov_group_idx = 2;
    fix_group_idx = 1;
    map_group_name_2_idx = containers.Map(image_group_list, {1, 2});
    mov_group = image_group_list{mov_group_idx};
    fix_group = image_group_list{fix_group_idx};
    %% Match ROI - assume roi index does not match
    image_group_info_cell = cell(num_image_group, 1);
    for iter_group = 1 : num_image_group
        image_group = image_group_list{iter_group};
        image_group_info_cell{iter_group} = DataManager.load_data(fullfile(...
            DataManager.fp_raw_data_folder(dataset_name, stack), image_group, ...
            sprintf('image_group_info.mat')));
    end
    assert(image_group_info_cell{1}.num_stack == image_group_info_cell{2}.num_stack);
    %% Load all the registration result in one stack
    num_roi = image_group_info_cell{1}.num_stack ;
    data_cell = cell(num_roi, 1);
    graph_cell = cell(2, num_roi);
    for iter_roi = 1 : num_roi
        registered_str = DataManager.load_data(fullfile(DataManager.fp_analysis_data_folder(dataset_name, ...
            stack), 'Matched_vessel', sprintf('%s_%s_%d_%s_%s_to_%s.mat', ...
                dataset_name, stack, iter_roi, 'rigid', ...
                mov_group, fix_group)));
        for iter_im_group = 1 : num_image_group
            tmp_im_group = image_group_list{iter_im_group};
            switch tmp_im_group
                case 'In_vivo'
                    tmp_grid_ver = sprintf('240_cube_%s_%d_re', tmp_im_group, iter_roi);
                case 'Post_perfusion'
                    tmp_grid_ver = sprintf('240_cube_%s_%d', tmp_im_group, iter_roi);
            end
            tmp_graph = DataManager.load_graph_in_block(dataset_name, stack, ...
                tmp_grid_ver, 0, 0, 0);
            %% Post processing
            tmp_graph.link.cc_com_sub = cellfun(@(x) mean(fun_ind2sub(...
                tmp_graph.num.mask_size, x), 1), tmp_graph.link.cc_ind, ...
                'UniformOutput', false);
            tmp_graph.link.cc_com_sub = cat(1, tmp_graph.link.cc_com_sub{:});
            tmp_graph.link.features.cc_com_sub_3 = tmp_graph.link.cc_com_sub(:, 3);
            tmp_graph.link.features.num_voxel = cellfun(@numel, tmp_graph.link.cc_ind);
            tmp_graph.link.features.radius_cv = tmp_graph.link.features.radius_std ./ ...
                tmp_graph.link.features.radius_mean;          
            
            switch tmp_im_group
                case 'In_vivo'
                    % Further select valid parameter estimations to refine
                    % radius estimation 
                    tmp_valid_radius_c = cellfun(@(x1, x2) x1(x2 > 0.8, :), ...
                        tmp_graph.link.cc_r, tmp_graph.link.cc_est_corr, ...
                        'UniformOutput', false);      
                    tmp_graph.link.features.num_valid_est = cellfun(@numel, tmp_valid_radius_c);
                    tmp_graph.link.features.valid_radius_median = cellfun(@(x) median(x, 1, 'omitnan'), ...
                        tmp_valid_radius_c, 'UniformOutput', false);
                    tmp_graph.link.features.valid_radius_median = cat(1, tmp_graph.link.features.valid_radius_median{:});
                    
                    tmp_graph.link.features.valid_radius_mean = cellfun(@(x) mean(x, 1, 'omitnan'), ...
                        tmp_valid_radius_c, 'UniformOutput', false);
                    tmp_graph.link.features.valid_radius_mean = cat(1, tmp_graph.link.features.valid_radius_mean{:});
                    
                    tmp_graph.link.features.valid_radius_std = cellfun(@(x) std(x, 0, 1, 'omitnan'), ...
                        tmp_valid_radius_c, 'UniformOutput', false);
                    tmp_graph.link.features.valid_radius_std = cat(1, tmp_graph.link.features.valid_radius_std{:});
                    
                    tmp_graph.link.features.valid_radius_cv = tmp_graph.link.features.valid_radius_std ./ ...
                        tmp_graph.link.features.valid_radius_mean;
                case 'Post_perfusion'
                    
            end            
            graph_cell{iter_im_group, iter_roi} = tmp_graph;
        end
        %% Add computation
        registered_str = fun_analysis_post_process_registered_vessel_graph(...
            registered_str, graph_cell{map_group_name_2_idx(registered_str.Fixed_image_group), iter_roi}, ...
            graph_cell{map_group_name_2_idx(registered_str.Moving_image_group), iter_roi});
             
        data_cell{iter_roi} = registered_str;
    end
    stack_data{iter_stack} = data_cell;
    fprintf('Finish loading data for stack %s\n', stack);
end
%% Post processing
vis_stack_group = {1, 2, 3, 4, 5, 6, 1:3, 1:4, 1:6};
num_vis_group = numel(vis_stack_group);
max_r = 10; % um
calibration_data = cell(num_vis_group, 1);
fixed_im_group = stack_data{1}{1}.Fixed_image_group;
moving_im_group =stack_data{1}{1}.Moving_image_group;

merge_stack_name = 'DK202004_merged';
visualization_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, merge_stack_name), ...
    'stack_statistics');
spline_skip_step = 2;
vis_stat_name = 'median';
equal_r_th = 8;
flat_r_th = 0.9;
min_datapont_fit = 20;
for iter_vis_group = 1 : num_vis_group
    tmp_calibration_str = struct;
    tmp_vis_stack_idx = vis_stack_group{iter_vis_group};
    num_vis_stack = numel(tmp_vis_stack_idx);
    merged_data = cat(1, stack_data{tmp_vis_stack_idx});
    merged_data = fun_analysis_radius_calibration_merge_registeration(merged_data);
    
    tmp_calibration_str.stack_list = stack_list(tmp_vis_stack_idx);
    tmp_calibration_str.Fixed_image_group = fixed_im_group;
    tmp_calibration_str.Moving_image_group = moving_im_group;
    %% Registration between vessel segments
    fixed_r = merged_data.cc.Fixed.radius_median;
    moving_r = merged_data.cc.Moving.radius_median;
    % Selection 
    n_diff_dt_N_prof_c = fun_normalized_difference(fixed_r(:, 2), fixed_r(:, 3), false);
    invivo_plot_idx = 2;
    selected_Q = merged_data.cc.Fixed.radius_cv(:, invivo_plot_idx) < 0.2 & merged_data.cc.Moving.radius_cv < 0.2 & ...
        merged_data.cc.Fixed.num_voxel > 10 & merged_data.cc.Moving.num_voxel > 10 & ...
        merged_data.cc.Fixed.cc_com_sub_3 < inf & merged_data.cc.Moving.cc_com_sub_3 < inf & ...
        abs(n_diff_dt_N_prof_c) < 0.2 & merged_data.cc.Dist_cc_mean < 4;

    tmp_str = fun_analysis_radius_calibration_str(fixed_r(selected_Q, invivo_plot_idx), ...
        moving_r(selected_Q), vis_stat_name, [0.5, 10], 25, flat_r_th, equal_r_th, 2, min_datapont_fit);
    tmp_str.Fixed_image_group = fixed_im_group;
    tmp_str.Moving_image_group = moving_im_group;
    % Save result
    tmp_calibration_str.cc = tmp_str;
    % Visualization
    [fig_hdl, ax_hdl] = fun_analysis_radius_calibration_vis(tmp_str, vis_stat_name);
    vis_x = linspace(tmp_str.Hist_edge_range(1), tmp_str.Hist_edge_range(2), 50);
    hold(ax_hdl, 'on');
    spline_hdl = plot(ax_hdl, vis_x, tmp_str.spline_itp(vis_x), 'LineWidth', 1.5);
    formula_hdl = plot(ax_hdl, vis_x, tmp_str.formula(vis_x), 'LineWidth', 1.5);
    
    leg_str = sprintf('r_i = r_p + r_0/(1 + r_p)\nr_0 = %.2f\nR^2-Adjusted = %.2f', ...
        tmp_str.formula_r0, tmp_str.formula_R2);
    leg_hdl = legend([formula_hdl, spline_hdl], leg_str, sprintf('Spline (%.2f, %.2f)', flat_r_th, equal_r_th), 'Location', 'southeast');
    
    if num_vis_stack == 1
        title_head_name = strrep(stack_list{tmp_vis_stack_idx}, '_', ' ');
    else
        title_head_name = sprintf('%d mice', num_vis_stack);
    end
    ax_hdl.Title.String = sprintf('%s, %d segments', title_head_name, tmp_str.Num_segments);
    %% Write figure
    fig_fp = fullfile(visualization_folder, sprintf('%s_%s_vessel_segment_%s_radius_comparison.png', ...
        dataset_name, title_head_name, vis_stat_name));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    delete(fig_hdl);
    %% Registration between skeleton voxels
    fixed_r = merged_data.voxel.Fixed_r_um;
    moving_r = merged_data.voxel.Moving_r_um;
    % Selection 
    voxel_r_n_diff = fun_normalized_difference(fixed_r(:, 2), fixed_r(:, 3));
    voxel_r_n_diff_dt_r0 = fun_normalized_difference(fixed_r(:, 1), fixed_r(:, 2));

    selected_Q = abs(voxel_r_n_diff) < 0.1 &...
        abs(voxel_r_n_diff_dt_r0 - median(voxel_r_n_diff_dt_r0, 'omitnan')) < 0.1 & ...
        merged_data.voxel.Fixed_est_corr > 0.95 & ...
        merged_data.voxel.Dist_fixed_2_moving < 4;
        
    tmp_str = fun_analysis_radius_calibration_str(fixed_r(selected_Q, invivo_plot_idx), ...
        moving_r(selected_Q), vis_stat_name, [0.5, 10], 30, flat_r_th, equal_r_th, 2, min_datapont_fit);
    tmp_str.Fixed_image_group = fixed_im_group;
    tmp_str.Moving_image_group = moving_im_group;
    % Save result
    tmp_calibration_str.voxel = tmp_str;
    % Visualization 
    [fig_hdl, ax_hdl] = fun_analysis_radius_calibration_vis(tmp_str, vis_stat_name);
    vis_x = linspace(tmp_str.Hist_edge_range(1), tmp_str.Hist_edge_range(2), 50);
    hold(ax_hdl, 'on');
    spline_hdl = plot(ax_hdl, vis_x, tmp_str.spline_itp(vis_x), 'LineWidth', 1.5);
    formula_hdl = plot(ax_hdl, vis_x, tmp_str.formula(vis_x), 'LineWidth', 1.5);
    
    leg_str = sprintf('r_i = r_p + r_0/(1 + r_p)\nr_0 = %.2f\nR^2-Adjusted = %.2f', ...
        tmp_str.formula_r0, tmp_str.formula_R2);
    leg_hdl = legend([formula_hdl, spline_hdl], leg_str, sprintf('Spline (%.2f, %.2f)', flat_r_th, equal_r_th), 'Location', 'southeast');
    
    % Visualization
    if num_vis_stack == 1
        title_head_name = strrep(stack_list{tmp_vis_stack_idx}, '_', ' ');
    else
        title_head_name = sprintf('%d mice', num_vis_stack);
    end
    ax_hdl.Title.String = sprintf('%s, %d voxles', title_head_name, tmp_str.Num_segments);
    %% Write figure
    fig_fp = fullfile(visualization_folder, sprintf('%s_%s_skeleton_voxel_%s_radius_comparison.png', ...
        dataset_name, title_head_name, vis_stat_name));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    delete(fig_hdl);
    %% Record result
    calibration_data{iter_vis_group} = tmp_calibration_str;
end
calibration_result_fp = fullfile(DataManager.fp_metadata_file(dataset_name, ...
    merge_stack_name, sprintf('%s_%s_%s_calibration_curve.mat', dataset_name, merge_stack_name, vis_stat_name)));
DataManager.write_data(calibration_result_fp, calibration_data);
%%
calibration_data = DataManager.load_data(fullfile(DataManager.fp_metadata_file(dataset_name, ...
    merge_stack_name, sprintf('%s_%s_%s_calibration_curve.mat', dataset_name, merge_stack_name, 'mean'))));
% Overlay 3 mice, 4 mice and 6 mice
vis_itp_x = linspace(0.5, 10, 50);
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
for iter_vis = 7 : 9
    tmp_str = calibration_data{iter_vis};
    errorbar(ax_hdl, tmp_str.voxel.Fixed_binned_by_Moving.x_bin_val, ...
        tmp_str.voxel.Fixed_binned_by_Moving.y_mean, ...
        tmp_str.voxel.Fixed_binned_by_Moving.y_std, 'LineWidth', 1.5);
    hold(ax_hdl, 'on');
end
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.XLim = [0.5, 10];
ax_hdl.YLim = [0.5, 10];
ax_hdl.XScale = 'log';
ax_hdl.YScale = 'log';
ax_hdl.XTick = [0.5 : 0.5 : 1.5, 2, 3, 4 : 2 : 10];
ax_hdl.YTick = [0.5 : 0.5 : 1.5, 2, 3, 4 : 2 : 10];
grid(ax_hdl, 'on');
plot(ax_hdl, 0 : 0.1 : 10, 0 : 0.1 : 10, 'k-.', 'LineWidth', 1);
ax_hdl.XLabel.String = sprintf('%s radius (\\mum)', strrep(tmp_str.Moving_image_group, '_', ' '));
ax_hdl.YLabel.String = sprintf('%s radius (\\mum)', strrep(tmp_str.Fixed_image_group, '_', ' '));
leg_hdl = legend('3 mice', '4 mice', '6 mice', 'Location', 'northwest');
ax_hdl.FontSize = 14;

fig_fp = fullfile(visualization_folder, sprintf('%s_%s_skeleton_voxel_%s_radius_calibration_curves.png', ...
    dataset_name, merge_stack_name, vis_stat_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
delete(fig_hdl);