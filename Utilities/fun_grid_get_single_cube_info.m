function cube_info = fun_grid_get_single_cube_info(grid_info, cube_label)

cube_info = struct;
% General information
cube_info.dataset_name = grid_info.dataset_name;
cube_info.stack = grid_info.stack;
cube_info.dataset_size = grid_info.data_size;
cube_info.voxel_size_um = grid_info.voxel_size_um;
cube_info.grid_version = grid_info.version;
cube_info.dataset_grid_size = grid_info.grid_size;
% Cube specific information
cube_info.grid_label = cube_label;
cube_info.global_bbox_mmxx = grid_info.bbox_xyz_mmxx_list(cube_label, :);
cube_info.global_bbox_mmll = grid_info.bbox_xyz_mmll_list(cube_label, :);
cube_info.idx_1 = grid_info.bbox_grid_sub_list(cube_label, 1);
cube_info.idx_2 = grid_info.bbox_grid_sub_list(cube_label, 2);
cube_info.layer = grid_info.bbox_grid_sub_list(cube_label, 3);
cube_info.grid_ind = grid_info.bbox_grid_ind_list(cube_label);
cube_info.block_size = cube_info.global_bbox_mmll(4:6);
end