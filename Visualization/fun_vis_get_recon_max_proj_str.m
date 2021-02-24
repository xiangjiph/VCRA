function proj_im_str = fun_vis_get_recon_max_proj_str(wb_max_proj_str, vis_proj_direction)

grid_info = wb_max_proj_str.grid_info;

vis_plan_dim = setdiff([1,2,3], vis_proj_direction);
vis_proj_plane_name_list = {'horizontal', 'sagittal', 'coronal'};
vis_proj_plane_name = vis_proj_plane_name_list{vis_proj_direction};
layer_list = grid_info.layer;
starting_idx_3 = layer_list(1);
end_idx_3 = layer_list(end);

vis_mask_cell = cell(wb_max_proj_str.grid_size);
switch vis_proj_direction
    case 1
        vis_perm_order = [2,3,1];        
        vis_mask_cell(wb_max_proj_str.grid_ind) = wb_max_proj_str.mask_max_proj_1;
        grid_ll_1 = squeeze(grid_info.grid2D.ll_array(2, 1, :));
        grid_ll_2 = grid_info.gridZ.ll(starting_idx_3 : end_idx_3);
        [tmp_bbox_size_list_1, tmp_bbox_size_list_2] = ndgrid(grid_ll_1, grid_ll_2);
        proj_plane_grid_ll_list = cat(2, tmp_bbox_size_list_1(:), tmp_bbox_size_list_2(:));
        vis_sub_min = [1, grid_info.gridZ.mmxx(starting_idx_3, 1)];
    case 2
        vis_perm_order = [1,3,2];
        vis_mask_cell(wb_max_proj_str.grid_ind) = wb_max_proj_str.mask_max_proj_2;
        grid_ll_1 = squeeze(grid_info.grid2D.ll_array(1, :, 1));
        grid_ll_2 = grid_info.gridZ.ll(starting_idx_3 : end_idx_3);
        [tmp_bbox_size_list_1, tmp_bbox_size_list_2] = ndgrid(grid_ll_1, grid_ll_2);
        proj_plane_grid_ll_list = cat(2, tmp_bbox_size_list_1(:), tmp_bbox_size_list_2(:));
        vis_sub_min = [1, grid_info.gridZ.mmxx(starting_idx_3, 1)];
    case 3
        vis_perm_order = [1,2,3];
        vis_mask_cell(wb_max_proj_str.grid_ind) = wb_max_proj_str.mask_max_proj_3;
        proj_plane_grid_ll_list = grid_info.grid2D.ll;
        vis_sub_min = [1, 1];
end
vis_mask_proj = permute(vis_mask_cell, vis_perm_order);
vis_subgrid_bbox_label_array = grid_info.bbox_grid_label_array(:, :, layer_list);
vis_subgrid_bbox_label_array = permute(vis_subgrid_bbox_label_array, vis_perm_order);
% Stitch max projection of cubes into one matrix for the entire layer
vis_im_cell = cell(grid_info.grid_size(vis_proj_direction), 1);
for iter_layer = 1 : grid_info.grid_size(vis_proj_direction)
    vis_im_cell{iter_layer} = fun_stitch_patchs_with_overlap_with_bbox_info(vis_mask_proj(:, :, iter_layer), ...
        proj_plane_grid_ll_list, [32, 32], 'logical');
end

proj_im_str.vis_im_cell = vis_im_cell;
proj_im_str.vis_proj_plane_name = vis_proj_plane_name;
proj_im_str.vis_sub_min = vis_sub_min;
proj_im_str.vis_subgrid_bbox_label_array = vis_subgrid_bbox_label_array;
proj_im_str.vis_perm_order = vis_perm_order;
proj_im_str.vis_plan_dim = vis_plan_dim;
proj_im_str.vis_proj_direction = vis_proj_direction;
end