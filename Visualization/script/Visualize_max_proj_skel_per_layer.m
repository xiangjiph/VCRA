set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'mouselight_1';
skl_grid_name = '240_cube_auto';
image_grid_name = '240_cube';
grid_c_version = '240_cube_combined_5';
grid_c_info = DataManager.load_grid(dataset_name, stack, grid_c_version);
grid_info = grid_c_info.grid_ori;
layer_list = 9 : grid_info.num_grid_layer;
num_layer_to_process = numel(layer_list);
internal_offset = 16;
im_size_xy = grid_info.data_xy_size;
%% Load skeleton and max project to the xy plane
pool_obj = gcp('nocreate');
delete(pool_obj);
parpool(8);
wb_skel_proj_3 = cell(num_layer_to_process, 1);
parfor iter_layer = 1 : num_layer_to_process
    tmp_tic = tic;
    grid_layer = layer_list(iter_layer);
    layer_valid_bbox_idx_mat = grid_info.bbox_xy_linear_idx_mat{grid_layer};
    block_grid_sub = grid_info.bbox_grid_sub{grid_layer};
    layer_valid_bbox_idx = grid_info.bbox_xy_idx{grid_layer};
    num_valid_graph = size(block_grid_sub, 1);
    layer_valid_bbox_xy_mmxx = grid_info.bbox_xy_mmxx{grid_layer};    
    layer_skel_proj_ind_cell = cell(num_valid_graph, 1);
    for iter_block = 1 : num_valid_graph
        tmp_idx_1 = block_grid_sub(iter_block, 1);
        tmp_idx_2 = block_grid_sub(iter_block, 2);
        tmp_layer = block_grid_sub(iter_block, 3);
        try
            tmp_skl_str = DataManager.load_block_skl(dataset_name, stack, skl_grid_name, ...
                tmp_idx_1, tmp_idx_2, tmp_layer);
            tmp_skl_sub = fun_ind2sub(tmp_skl_str.block_size, tmp_skl_str.ind);
            tmp_skl_sub_global = bsxfun(@plus, tmp_skl_sub, tmp_skl_str.global_bbox_mmll(1:3) - 1);
            % Indices of the skeleton voxels projected to xy plane
            tmp_skl_ind_global_2d = sub2ind(im_size_xy, tmp_skl_sub_global(:, 1), tmp_skl_sub_global(:, 2));
            layer_skel_proj_ind_cell{iter_block} = tmp_skl_ind_global_2d;
        catch ME
            fprintf('Unable to load reconstructed mask (%d, %d, %d)\n', tmp_idx_1, tmp_idx_2, tmp_layer);
            continue;
        end
    end
    tmp_sec_skel_mask = false(im_size_xy);
    tmp_sec_skel_mask(cat(1, layer_skel_proj_ind_cell{:})) = true;
    wb_skel_proj_3{iter_layer} = tmp_sec_skel_mask;
    fprintf('Finish processing layer %d. Elapsed time is %f seconds\n', grid_layer, toc(tmp_tic));
end
%%
wb_skel_proj_3_stack = cat(3, wb_skel_proj_3{:});
implay(wb_skel_proj_3_stack);
DataManager.write_tiff_stack(wb_skel_proj_3_stack(:, :, 29),...
    DataManager.fp_constructor(dataset_name, stack, 'visualization', 'WholeBrain_mouselight_1_240_cube_max_projection_layer_29.tiff'));

% avi_str = VideoWriter(DataManager.fp_constructor(dataset_name, stack, 'visualization', 'WholeBrain_mouselight_1_240_cube_max_projection.avi',true));
% avi_str.FrameRate = 5;
% avi_str.Quality = 75;
% open(avi_str);
% for frame_idx = 1 : size(wb_skel_proj_3_stack, 3)
%     writeVideo(avi_str, uint8(wb_skel_proj_3_stack(:,:,frame_idx)) * 255);
% end
% close(avi_str)