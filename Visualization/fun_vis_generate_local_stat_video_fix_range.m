function exit_code = fun_vis_generate_local_stat_video_fix_range(data_array, plot_name, colorbar_name, grid_info, proj_im_str, cbar_limit_low, cbar_limit_high, save_im_Q)

if nargin < 8
    save_im_Q = false;
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
vis_data_array(grid_info.bbox_grid_ind_list) = data_array;
vis_data_array = permute(vis_data_array, vis_perm_order);

vis_opt = struct;

vis_opt.folder_name = plot_name;
vis_opt.fig_name = plot_name;
vis_opt.color_bar_label = colorbar_name;

vis_opt.pixel_size = vis_downsample_rate;
vis_opt.cbar_low = cbar_limit_low;
vis_opt.cbar_high = cbar_limit_high;
vis_opt.vis_figQ = false;
vis_opt.return_frameQ = true;
vis_opt.cbar_tick_lable = cellfun(@(x) num2str(x, '%.2e'),...
    num2cell(linspace(vis_opt.cbar_low, vis_opt.cbar_high, 11)), 'UniformOutput', false);
vis_opt.video_file_name = sprintf('%s_%s_%s_%s.avi', dataset_name, stack, grid_info.version, vis_opt.fig_name);
vis_opt.video_file_path = DataManager.fp_visualization_image(dataset_name, stack, vis_opt.video_file_name, fullfile(vis_folder_name, vis_opt.folder_name));

avi_str = VideoWriter(vis_opt.video_file_path);
avi_str.FrameRate = 5;
avi_str.Quality = 80;
num_layer_to_plot = size(vis_data_array, 3);
open(avi_str);
for layer_idx = 1 : num_layer_to_plot
    fprintf('Writing frame %d/%d\n', layer_idx, num_layer_to_plot);
    vis_data_mat = vis_data_array(:, :, layer_idx);
    % Color bar saturation level:
    
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
    
    patch_value = vis_data_mat(vis_layer_patch_idx);
    layer_frame = fun_vis_image_with_patch(vis_im, patch_value, vis_patch_bbox_xy, vis_opt);
    writeVideo(avi_str, layer_frame);
end
close(avi_str)
exit_code = 0;
end