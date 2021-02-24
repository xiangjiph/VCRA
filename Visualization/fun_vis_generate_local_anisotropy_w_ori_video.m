function exit_code = fun_vis_generate_local_anisotropy_w_ori_video(data_list, ...
    plot_name, color_bar_name, ori_vec_list,  grid_info, proj_im_str, ...
    vis_method, save_im_Q, ptl_l, ptl_h)

if nargin < 7
    save_im_Q = false;
    ptl_l = 10;
    ptl_h = 90;
elseif nargin < 8
    ptl_l = 10;
    ptl_h = 90;
end

vis_downsample_rate = 4;

persistent DataManager
if isempty(DataManager)
    DataManager = FileManager;
end

dataset_name = grid_info.dataset_name;
stack = grid_info.stack;
grid_version = grid_info.version;
vis_perm_order = proj_im_str.vis_perm_order;
vis_subgrid_bbox_label_array = proj_im_str.vis_subgrid_bbox_label_array;
vis_im_cell = proj_im_str.vis_im_cell;
vis_plan_dim = proj_im_str.vis_plan_dim;
vis_sub_min = proj_im_str.vis_sub_min;

vis_folder_name = sprintf('%s_local_stat_%s', grid_version, proj_im_str.vis_proj_plane_name);
% Construct array and permute the array. 
vis_data_array = nan(grid_info.grid_size);
vis_data_array(grid_info.bbox_grid_ind_list) = data_list;
vis_data_array = permute(vis_data_array, vis_perm_order);

ori_vec = nan([3, grid_info.grid_size]);
ori_vec(:, grid_info.bbox_grid_ind_list) = ori_vec_list';
ori_vec_cell = {squeeze(ori_vec(1, :, :, :)), squeeze(ori_vec(2, :, :, :)), squeeze(ori_vec(3, :, :, :))};
ori_vec_cell = cellfun(@(x) permute(x, vis_perm_order), ori_vec_cell, 'UniformOutput', false);
ori_vec_vis_1 = ori_vec_cell{vis_plan_dim(1)};
ori_vec_vis_2 = ori_vec_cell{vis_plan_dim(2)};

% Fix the color bar limit for the whole brain
vis_valid_array_Q = ~isnan(vis_data_array);

vis_data_valid = vis_data_array(vis_valid_array_Q);
ori_vec_vis_1(~vis_valid_array_Q) = nan;
ori_vec_vis_2(~vis_valid_array_Q) = nan;

cbar_limit_low = prctile(vis_data_valid, ptl_l);
cbar_limit_high = prctile(vis_data_valid, ptl_h);
vis_opt = struct;

vis_opt.folder_name = sprintf('%s_%s', plot_name, vis_method);
vis_opt.fig_name = sprintf('%s_%s', plot_name, vis_method);
vis_opt.color_bar_label = color_bar_name;

vis_opt.pixel_size = vis_downsample_rate;
vis_opt.cbar_low = cbar_limit_low;
vis_opt.cbar_high = cbar_limit_high;
vis_opt.vis_figQ = false;
vis_opt.return_frameQ = true;
vis_opt.cbar_tick_lable = cellfun(@(x) num2str(x, '%.1e'),...
    num2cell(linspace(vis_opt.cbar_low, vis_opt.cbar_high, 11)), 'UniformOutput', false);
vis_opt.video_file_name = sprintf('%s_%s_%s_%s.avi', dataset_name, stack, grid_info.version, vis_opt.fig_name);
vis_opt.video_file_path = DataManager.fp_visualization_image(dataset_name, stack, vis_opt.video_file_name, fullfile(vis_folder_name, vis_opt.folder_name));

switch vis_method
    case 'vec_proj'
        vis_opt.cmap_name = 'winter';
    case 'color_coded'
        vis_opt.cmap_name = 'hue';
end

avi_str = VideoWriter(vis_opt.video_file_path);
avi_str.FrameRate = 5;
avi_str.Quality = 90;
num_layer_to_plot = size(vis_data_array, 3);
open(avi_str);
for layer_idx = 1 : num_layer_to_plot
    fprintf('Writing frame %d/%d\n', layer_idx, num_layer_to_plot);
    vis_data_mat = vis_data_array(:, :, layer_idx);    
    vis_opt.fig_title = sprintf('%s in grid layer %d', strrep(vis_opt.fig_name, '_', ' '), layer_idx);
    vis_opt.file_name = sprintf('%s_%s_%s_%s_layer_%d.png', dataset_name, stack, grid_info.version, vis_opt.fig_name, layer_idx);
    if save_im_Q
        vis_opt.file_path = DataManager.fp_visualization_image(dataset_name, stack, vis_opt.file_name, fullfile(vis_folder_name, vis_opt.folder_name));
    end
    % Determine the position of the patches
    vis_layer_patch_idx = find(~isnan(vis_data_mat));
    vis_layer_patch_label_mat = vis_subgrid_bbox_label_array(:, :, layer_idx);
    vis_layer_patch_label = vis_layer_patch_label_mat(vis_layer_patch_idx);
    vis_layer_bbox_mmxx = grid_info.bbox_xyz_mmxx_list(vis_layer_patch_label, :);
    
    vis_im = vis_im_cell{layer_idx};
    vis_im_size = ceil(size(vis_im)./ vis_opt.pixel_size);
    vis_im = imresize(vis_im, vis_im_size, 'Method', 'nearest');
    vis_patch_bbox_xy = vis_layer_bbox_mmxx(:, [vis_plan_dim, vis_plan_dim + 3]);
    vis_patch_bbox_xy = vis_patch_bbox_xy - [vis_sub_min, vis_sub_min] + 1;
    % Scale the patches
    vis_patch_bbox_xy = ceil(vis_patch_bbox_xy ./ vis_opt.pixel_size);
    
    switch vis_method
        case 'vec_proj'
            % Compute the 2D orientation vector
            tmp_ori_vec_2d_1 = ori_vec_vis_1(:, :, layer_idx);
            tmp_ori_vec_2d_2 = ori_vec_vis_2(:, :, layer_idx);
            tmp_ori_vec_1 = tmp_ori_vec_2d_1(vis_layer_patch_idx);
            tmp_ori_vec_2 = tmp_ori_vec_2d_2(vis_layer_patch_idx);
            ori_vec_2d_list = cat(2, tmp_ori_vec_1, tmp_ori_vec_2);
            
            patch_value = vis_data_mat(vis_layer_patch_idx);
            layer_frame = fun_vis_single_section_local_network_orientation(vis_im, patch_value, ori_vec_2d_list, vis_patch_bbox_xy, vis_opt);
        case 'color_coded'
            tmp_ori_vec_3d_1 = ori_vec_cell{1}(:, :, layer_idx);
            tmp_ori_vec_3d_2 = ori_vec_cell{2}(:, :, layer_idx);
            tmp_ori_vec_3d_3 = ori_vec_cell{3}(:, :, layer_idx);
            tmp_ori_vec_3d_1 = tmp_ori_vec_3d_1(vis_layer_patch_idx);
            tmp_ori_vec_3d_2 = tmp_ori_vec_3d_2(vis_layer_patch_idx);
            tmp_ori_vec_3d_3 = tmp_ori_vec_3d_3(vis_layer_patch_idx);
            
            ori_vec_3d_list = cat(2, tmp_ori_vec_3d_3, tmp_ori_vec_3d_2, tmp_ori_vec_3d_1);            
            layer_frame = fun_vis_single_section_local_network_orientation_by_color_code(vis_im, ori_vec_3d_list, vis_patch_bbox_xy, vis_opt);
    end
    writeVideo(avi_str, layer_frame);
end
close(avi_str)
exit_code = 0;
end