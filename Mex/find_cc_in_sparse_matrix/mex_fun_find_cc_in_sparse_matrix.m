function [PixelIdxList, num_cc] = mex_fun_find_cc_in_sparse_matrix(idx_list, neighbor_voxel_idx)

num_voxel = length(idx_list);
cc_length = 10000;
num_cc = 0;
PixelIdxList = cell(1,num_voxel);

% Initialization
t_node_unvisited = true(num_voxel,1);
t_start_search_idx = 1;
t_num_unvisited_points = num_voxel;
queue = zeros(cc_length,1);
t_cc_node_ind_list = zeros(cc_length,1);
t_current_id = 0;
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
            t_neigh_idx = neighbor_voxel_idx(:, t_current_id);                       
            t_neigh_idx = t_neigh_idx(t_neigh_idx>0);
            t_neigh_idx = t_neigh_idx(t_node_unvisited(t_neigh_idx));
            t_num_effect_neighbor = length(t_neigh_idx);
            
            if t_num_effect_neighbor == 0
                continue;
            else
                % Add the new voxels index into the queue and the cc list
                t_node_unvisited(t_neigh_idx) = false;
                t_cc_node_ind_list(t_num_vol+1 : t_num_vol + t_num_effect_neighbor) = idx_list(t_neigh_idx);
                queue(queue_end_pointer+1 : queue_end_pointer + t_num_effect_neighbor) = t_neigh_idx;
                queue_end_pointer = queue_end_pointer + t_num_effect_neighbor;
                t_num_vol = t_num_vol + t_num_effect_neighbor;
            end
        end
        num_cc = num_cc + 1;
        PixelIdxList{num_cc} = t_cc_node_ind_list(1:t_num_vol);
        t_num_unvisited_points = t_num_unvisited_points - t_num_vol;
    end
end
end
