function cc = fun_cc_in_sparse_matrix_v1(idx_list, mask_size)
% This function collect the connected components in a mask 
% The mask is specified by two inputs:
%   idx_list: N-by-1 double precision array. List of mask voxel indices
%   mask_size: size of the mask
% Output:
%   cc: cell array, each contains the indics in the idx_list
% 
% For small matrix, this function performs less preferable than coverting
% the sparse matrix to the full matric and use the built-in bwconcomp.
% Need to test the performance on much larger array later.
% This function is needed for finding the connected components when the
% mask array is very large and very sparse that the full array cannot fit
% into the memory. 

% Written by Xiang Ji on Aug 23, 2018

% Generate 26 neighbor indices addition coefficient
mask_size_pad = mask_size + 2;
[tmp1, tmp2, tmp3] = ndgrid(1:3);
neighbor_add_coeff = sub2ind(mask_size_pad, tmp1(:), tmp2(:), tmp3(:));
neighbor_add_coeff = neighbor_add_coeff - neighbor_add_coeff(14);
neighbor_add_coeff(14) = [];
%% Pad array

[pos_1, pos_2, pos_3] = ind2sub(mask_size, idx_list);
idx_list_padded = sub2ind(mask_size_pad, pos_1 + 1, pos_2 + 1, pos_3 + 1); 
num_voxel = length(idx_list);
num_block_voxel_pad = prod(mask_size_pad);

sparse_matrix = sparse(idx_list_padded, ones(num_voxel, 1), ...
    1:num_voxel, num_block_voxel_pad, 1);
clear pos_1 pos_2 pos_3


% Internal parameters
    % Estimated maximum size of the connected components and the queue
cc_length = 1000;
t_num_cc = 0;
cc.PixelIdxList = cell(1,num_voxel);
t_node_unvisited = true(num_voxel,1);
t_start_search_idx = 1;
t_num_unvisited_points = num_voxel;
queue = zeros(cc_length,1);
t_cc_node_ind_list = zeros(cc_length,1);
while t_num_unvisited_points > 0
    for t_node_idx = t_start_search_idx : num_voxel
        if t_node_unvisited(t_node_idx)
            t_current_id = t_node_idx;
            t_start_search_idx = t_node_idx + 1;
            t_node_unvisited(t_node_idx) = false;
            break;
        end
    end
    
    if t_num_unvisited_points <= 0
        break;        
    else
        queue_start_pointer = 1;
        queue_end_pointer = 1;
        queue(queue_start_pointer) = t_current_id;
        t_cc_node_ind_list(1) = idx_list(t_current_id);
        t_num_vol = 1;
        while (queue_end_pointer - queue_start_pointer) >= 0
            t_current_id = queue(queue_start_pointer);
            queue_start_pointer = queue_start_pointer + 1;
            t_current_pos = idx_list_padded(t_current_id);            
            % Look for the 26 neighbors of the current point
            t_neigh_node_idx = full(sparse_matrix( t_current_pos + neighbor_add_coeff));
            t_neigh_node_idx = t_neigh_node_idx(t_neigh_node_idx>0);
            t_neigh_node_idx = t_neigh_node_idx(t_node_unvisited(t_neigh_node_idx));
            t_num_effect_neighbor = length(t_neigh_node_idx);
            
            if t_num_effect_neighbor == 0
                continue;
            else
                % Add the new voxels index into the queue and the cc list
                t_node_unvisited(t_neigh_node_idx) = false;
                t_cc_node_ind_list(t_num_vol+1 : t_num_vol + t_num_effect_neighbor) = idx_list(t_neigh_node_idx);
                queue(queue_end_pointer+1 : queue_end_pointer + t_num_effect_neighbor) = t_neigh_node_idx;
                queue_end_pointer = queue_end_pointer + t_num_effect_neighbor;
                t_num_vol = t_num_vol + t_num_effect_neighbor;
            end
        end
        t_num_cc = t_num_cc + 1;
        cc.PixelIdxList{t_num_cc} = t_cc_node_ind_list(1:t_num_vol);
        t_num_unvisited_points = t_num_unvisited_points - t_num_vol;
    end
end
cc.PixelIdxList = cc.PixelIdxList(1:t_num_cc);
cc.NumObjects = t_num_cc;
cc.Connectivity = 26;
cc.ImageSize = mask_size;
end
