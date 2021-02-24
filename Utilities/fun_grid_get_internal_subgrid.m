function grid_c = fun_grid_get_internal_subgrid(grid_c)
% This function is a helper function for generating the unbaised block
% statistics structure fron the reconstructed vessel graph. It adds several
% fields to the input combined grid structure, which will be used for
% determine from which combined grid the local statistics in each subgrid (
% 240 cube) should be computed. 
warning('This function need to be check, since sub_grid_ind is renamed to be sub_grid_label');

num_grid_c = size(grid_c.bbox_xyz_mmxx_grid_list, 1);
subgrid_size = size(grid_c.sub_grid_ind_array);
dist_to_boundary_cell = cell(subgrid_size);
% cube_valid_in_grid_c_label = zeros(num_cube, 1);
subgrid_valid_in_grid_c_label_array = nan(subgrid_size);
% disp('Find the combined grid label where the subgrid is an internal subgrid');
%
for iter_grid_c = 1 : num_grid_c
    grid_c_sub = grid_c.bbox_grid_sub_list(iter_grid_c, :);
    grid_c_bbox_xyz_mmll_grid = grid_c.bbox_xyz_mmll_grid_list(iter_grid_c, :);
%     tmp_subgrid_local_ind = grid_c.sub_grid_valid_idx{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)};
    tmp_subgrid_global_sub = grid_c.sub_grid_sub{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)};
    tmp_num_subgrid = size(tmp_subgrid_global_sub, 1);
        
    grid_c_bbox_size = grid_c_bbox_xyz_mmll_grid(4:6);
%     grid_c_local_valid_Q_array = false(grid_c_bbox_size);
%     grid_c_local_valid_Q_array(tmp_subgrid_local_ind) = true;
    % Compute distance to the boundary
%     dist_to_boudanry_pad = ~grid_c_local_valid_Q_array;
    dist_to_boudanry_pad = false(grid_c_bbox_size);
    dist_to_boudanry_pad = padarray(dist_to_boudanry_pad, [1,1,1], true, 'both');
    dist_to_boudanry_pad = bwdist(dist_to_boudanry_pad);
    grid_c_bbox_size_grid_pad = grid_c_bbox_size + 2;
    n26_ind_add_pad = fun_skeleton_neighbor_add_coeff_3D(grid_c_bbox_size_grid_pad, 26, false);
    
    tmp_subgrid_local_sub = tmp_subgrid_global_sub - grid_c_bbox_xyz_mmll_grid(1:3) + 1;
    tmp_subgrid_local_ind_pad = sub2ind(grid_c_bbox_size_grid_pad, tmp_subgrid_local_sub(:, 1) + 1, ...
        tmp_subgrid_local_sub(:, 2) + 1, tmp_subgrid_local_sub(:, 3) + 1);
    tmp_dist_to_boundary_27 = dist_to_boudanry_pad(bsxfun(@plus, tmp_subgrid_local_ind_pad', n26_ind_add_pad));
    % Truncate the distance to maximum 1.
    tmp_dist_to_boundary_27 = min(1, tmp_dist_to_boundary_27);
    for iter_subgrid = 1 : tmp_num_subgrid
        tmp_ind = sub2ind(subgrid_size, tmp_subgrid_global_sub(iter_subgrid, 1), ...
            tmp_subgrid_global_sub(iter_subgrid, 2), tmp_subgrid_global_sub(iter_subgrid, 3));
        tmp_dist = tmp_dist_to_boundary_27(:, iter_subgrid);
        if isempty(dist_to_boundary_cell{tmp_ind})
            dist_to_boundary_cell{tmp_ind} = tmp_dist;
            subgrid_valid_in_grid_c_label_array(tmp_ind) = iter_grid_c;
        else
            prev_dist_to_boundary = dist_to_boundary_cell{tmp_ind};
            if all(prev_dist_to_boundary <= tmp_dist) || sum(prev_dist_to_boundary) < sum(tmp_dist)
%                 disp('New label is better. Update. ');
                dist_to_boundary_cell{tmp_ind} = tmp_dist;
                subgrid_valid_in_grid_c_label_array(tmp_ind) = iter_grid_c;
%             elseif all(prev_dist_to_boundary >= tmp_dist)
%                 disp('Existing label is better');               
            else
%                 disp('Half-half. Do not update');
            end
        end
    end
end
%% 
grid_c.internal_subgrid_label_array = subgrid_valid_in_grid_c_label_array;
[grid_c.internal_subgrid_valid_Q, grid_c.internal_subgrid_valid_idx, ...
    grid_c.internal_subgrid_ind, grid_c.internal_subgrid_sub, grid_c.internal_subgrid_label,...
    grid_c.internal_subgrid_bbox_mmll, grid_c.internal_subgrid_bbox_mmxx] = deal(cell(grid_c.grid_size));
for iter_grid_c  = 1 : num_grid_c
    grid_c_sub = grid_c.bbox_grid_sub_list(iter_grid_c, :);
    tmp_subgrid_global_sub = grid_c.sub_grid_sub{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)};
    tmp_subgrid_global_ind = sub2ind(subgrid_size, tmp_subgrid_global_sub(:, 1), ...
        tmp_subgrid_global_sub(:, 2), tmp_subgrid_global_sub(:, 3));
    tmp_subgrid_valid_grid_label = subgrid_valid_in_grid_c_label_array(tmp_subgrid_global_ind);
    tmp_subgrid_is_validQ = (tmp_subgrid_valid_grid_label == iter_grid_c);
    tmp_valid_idx = grid_c.sub_grid_valid_idx{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)};
    
    grid_c.internal_subgrid_valid_Q{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)} = tmp_subgrid_is_validQ;
    grid_c.internal_subgrid_valid_idx{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)} =  tmp_valid_idx(tmp_subgrid_is_validQ);
    
    grid_c.internal_subgrid_label{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)} = ...
        grid_c.sub_grid_label{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)}(tmp_subgrid_is_validQ, :);
    
    grid_c.internal_subgrid_ind{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)} = ...
        grid_c.sub_grid_ind{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)}(tmp_subgrid_is_validQ, :);
    
    grid_c.internal_subgrid_sub{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)} = ...
        grid_c.sub_grid_sub{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)}(tmp_subgrid_is_validQ, :);
    
    grid_c.internal_subgrid_bbox_mmll{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)} = ...
        grid_c.sub_grid_bbox_mmll{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)}(tmp_subgrid_is_validQ, :);
    
    grid_c.internal_subgrid_bbox_mmxx{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)} = ...
        grid_c.sub_grid_bbox_mmxx{grid_c_sub(1), grid_c_sub(2), grid_c_sub(3)}(tmp_subgrid_is_validQ, :);
end
end