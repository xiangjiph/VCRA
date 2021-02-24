function [grid_c_skl_ind, grid_c_r, int_bbox_local_bbox_mmxx, int_bbox_exp_size]...
    = fun_simulation_OT_load_skeleton_in_grid(grid_c_info, load_skl_name, ...
    tmp_grid_c_ind, int_expand_half_length)

persistent DataManager
if isempty(DataManager)
    DataManager = FileManager([], true);
end
%%
if isempty(grid_c_info.internal_subgrid_valid_idx{tmp_grid_c_ind})
    fprintf('No valid internal sub-grid in this combined grid\n');
    [grid_c_skl_ind, grid_c_r, int_bbox_local_bbox_mmxx] = deal([]);
else
    %% Load all the subgrid data
    % All the subgrid inside the combined grid
    grid_c_subgrid_sub = grid_c_info.sub_grid_sub{tmp_grid_c_ind};
    
    % Internal subgrids inside the combined grid
    int_subgrid_bbox_mmxx_list = grid_c_info.internal_subgrid_bbox_mmxx{tmp_grid_c_ind};
    int_bbox_exp_mm = max([1,1,1], min(int_subgrid_bbox_mmxx_list(:, 1:3), [], 1) - int_expand_half_length);
    int_bbox_exp_xx = min(grid_c_info.data_size, max(int_subgrid_bbox_mmxx_list(:, 4:6), [], 1) + int_expand_half_length);
    int_bbox_exp_size = int_bbox_exp_xx - int_bbox_exp_mm + 1;

    % Load skeleton and combined them
    grid_c_num_subgrid = size(grid_c_subgrid_sub, 1);
    assert(nnz(grid_c_info.sub_grid_label_array{tmp_grid_c_ind} > 0) == grid_c_num_subgrid, 'Number of valid subgrid does not match');
    grid_c_skl_ind = cell(grid_c_num_subgrid, 1);
    grid_c_r = cell(grid_c_num_subgrid, 1);
    for iter_cube = 1 : grid_c_num_subgrid
        tmp_idx_1 = grid_c_subgrid_sub(iter_cube, 1);
        tmp_idx_2 = grid_c_subgrid_sub(iter_cube, 2);
        tmp_layer = grid_c_subgrid_sub(iter_cube, 3);
        try
            tmp_skel = DataManager.load_block_skl(grid_c_info.dataset_name,...
                grid_c_info.stack, load_skl_name, tmp_idx_1, tmp_idx_2, tmp_layer);
        catch
            fprintf('Fail to read the skeleton file of cube (%d, %d, %d) . Skip this block\n', tmp_idx_1, tmp_idx_2, tmp_layer);
            continue;
        end
        tmp_skl_sub = fun_ind2sub(tmp_skel.block_size, tmp_skel.ind);
        tmp_skl_sub_g = bsxfun(@plus, tmp_skl_sub, tmp_skel.global_bbox_mmll(1:3) - 1);
        tmp_skl_sub_int_exp = tmp_skl_sub_g - int_bbox_exp_mm + 1;
        
        tmp_in_bbox_Q = all(tmp_skl_sub_int_exp > 0, 2) & ...
            all( bsxfun(@le, tmp_skl_sub_int_exp, int_bbox_exp_size), 2);
        tmp_skl_sub_valid = tmp_skl_sub_int_exp(tmp_in_bbox_Q , : );
        
        grid_c_skl_ind{iter_cube} = sub2ind(int_bbox_exp_size, tmp_skl_sub_valid(:,1), ...
            tmp_skl_sub_valid(:, 2), tmp_skl_sub_valid(:,3));
        grid_c_r{iter_cube} = tmp_skel.r(tmp_in_bbox_Q);
    end
    grid_c_skl_ind = cat(1, grid_c_skl_ind{:});
    grid_c_r = cat(1, grid_c_r{:});
    assert(numel(grid_c_skl_ind) == numel(grid_c_r), 'Number of skeleton voxels does not equal to the number of radius');
    % Output
    [grid_c_skl_ind, tmp_unique_idx, ~] = unique(grid_c_skl_ind);
    grid_c_r = grid_c_r(tmp_unique_idx);
    int_bbox_local_bbox_mmxx = int_subgrid_bbox_mmxx_list - [int_bbox_exp_mm, int_bbox_exp_mm] + 1;
end
end