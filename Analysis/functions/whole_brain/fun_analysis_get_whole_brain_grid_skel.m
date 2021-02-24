function wb_skel_str = fun_analysis_get_whole_brain_grid_skel(dataset_name, stack, grid_info, skl_version)

persistent DataManager
if isempty(DataManager)
    DataManager = FileManager;
end
skel_ind_cell = cell(grid_info.grid_size);
skel_r_cell = cell(grid_info.grid_size);
pool_obj = gcp('nocreate');
delete(pool_obj);
parpool(16);
tmp_tic_1 = tic;
layer_list = grid_info.layer;
%% Parallel load skeleton
parfor iter_layer = 1 : numel(layer_list)
    tmp_tic = tic;
    grid_layer = layer_list(iter_layer);
    block_grid_sub = grid_info.bbox_grid_sub{grid_layer};
    num_valid_graph = grid_info.num_bbox_xy(grid_layer);
    tmp_ind_cell = cell(grid_info.grid_size(1:2));
    tmp_r_cell = cell(grid_info.grid_size(1:2));
    for iter_block = 1 : num_valid_graph
        tmp_idx_1 = block_grid_sub(iter_block, 1);
        tmp_idx_2 = block_grid_sub(iter_block, 2);
        tmp_layer = block_grid_sub(iter_block, 3);
        try
            tmp_skl = DataManager.load_block_skl(dataset_name, stack, skl_version, tmp_idx_1, tmp_idx_2, tmp_layer);
            tmp_ind = tmp_skl.ind;
            tmp_sub = fun_ind2sub(tmp_skl.block_size, tmp_ind);
            % Determine the valid bounding box
            tmp_bbox_mmxx = fun_get_block_valid_bbox_mmxx([tmp_skl.idx_1, tmp_skl.idx_2, tmp_skl.layer], ...
                grid_info.bbox_grid_label_array > 0, tmp_skl.block_size, grid_info.block_overlap);
            % Select the valid voxels in the bounding box            
            tmp_in_bbox_Q = fun_voxel_sub_in_bbox_mmxx_Q(tmp_sub, tmp_bbox_mmxx);
            tmp_sub = tmp_sub(tmp_in_bbox_Q, :);

            tmp_sub_global = bsxfun(@plus, tmp_sub, tmp_skl.global_bbox_mmll(1:3) - 1);
            tmp_ind_global = sub2ind(tmp_skl.dataset_size, tmp_sub_global(:, 1), tmp_sub_global(:, 2), tmp_sub_global(:, 3));
            
            tmp_ind_cell{tmp_idx_1, tmp_idx_2} = tmp_ind_global;
            tmp_r_cell{tmp_idx_1, tmp_idx_2} = tmp_skl.r(tmp_in_bbox_Q);
        catch
            fprintf('Failed to load skeleton in block (%d, %d, %d)\n', tmp_idx_1, tmp_idx_2, tmp_layer); 
        end
    end
    skel_ind_cell(:, :, iter_layer) = tmp_ind_cell;
    skel_r_cell(:, :, iter_layer) = tmp_r_cell;
    fprintf('Finish loading layer %d. Elapsed time is %f seconds\n', grid_layer, toc(tmp_tic));
end
pool_obj = gcp('nocreate');
delete(pool_obj);
%%
fprintf('Finish loading all the skeleton data. Elapsed time is %f seconds\n', toc(tmp_tic_1));
wb_skel_str.ind = skel_ind_cell;
wb_skel_str.r = skel_r_cell;
end