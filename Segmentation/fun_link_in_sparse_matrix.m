function cc = fun_link_in_sparse_matrix(link_idx_list, mask_size) %#codegen
% fun_link_in_sparse_matrix finds the connected component in sparse matrix.
% This algorithm is faster than fun_cc_in_sparse_matrix since it is
% designed specifically for skeleton segment tracing.
% Input:
%   idx_list: N-by-1 double precision array. List of mask voxel indices
%   mask_size: size of the mask
% Output:
%   cc: cell array, each contains the indics in the idx_list, which start
%   from one of the end of the link and end with the other end of the link
%   segment, which is not the case in fun_cc_in_sparse_matrix or bwconcomp

% Generate 26 neighbor indices addition coefficient
mask_size_pad = mask_size + 2;
[pos_1, pos_2, pos_3] = ndgrid(1:3);
neighbor_add_coeff = sub2ind(mask_size_pad, pos_1(:), pos_2(:), pos_3(:));
neighbor_add_coeff = neighbor_add_coeff - neighbor_add_coeff(14);
neighbor_add_coeff(14) = [];
% Pad array
[pos_1, pos_2, pos_3] = ind2sub(mask_size, link_idx_list);
idx_list_padded = sub2ind(mask_size_pad, pos_1 + 1, pos_2 + 1, pos_3 + 1); 
% clear pos_1 pos_2 pos_3
num_voxel = length(link_idx_list);
num_block_voxel_pad = prod(mask_size_pad);
% Find the position of the 26 neighbors for each node
if num_block_voxel_pad > 500^3 % This number should be optimized later. 
    % For large array, use sparse matrix
    mask_padded = sparse(idx_list_padded, ones(num_voxel, 1), ...
        1:num_voxel, num_block_voxel_pad, 1);
    neighbor_voxel_idx = full(mask_padded(bsxfun(@plus, neighbor_add_coeff, idx_list_padded')));
else
    mask_padded = zeros(num_block_voxel_pad, 1);
    mask_padded(idx_list_padded) = 1:num_voxel;
    neighbor_voxel_idx = mask_padded(bsxfun(@plus, neighbor_add_coeff, idx_list_padded'));
end
neighbor_voxel_idx = sort(neighbor_voxel_idx, 1, 'descend');
% Not sure if this will help. 
% neighbor_voxel_idx = neighbor_voxel_idx(1:2, :);

l_num_neighbor = sum(neighbor_voxel_idx>0, 1);
l_link_start_voxel_idx = find(l_num_neighbor==1);
l_iso_endpoint_voxel_idx = find(l_num_neighbor==0);
t_num_link_start_voxel = numel(l_link_start_voxel_idx);

t_link_idx = zeros(1000,1);
t_start_search_idx = 1;
% Single voxel link
cc.PixelIdxList = cell(1,num_voxel);
num_iso_endpoint = length(l_iso_endpoint_voxel_idx);
cc.PixelIdxList(1:num_iso_endpoint) = mat2cell(link_idx_list(l_iso_endpoint_voxel_idx), ones(num_iso_endpoint,1));
voxel_unvisited = true(num_voxel, 1);
voxel_unvisited(l_iso_endpoint_voxel_idx) = false;
% Link of length more than 1
num_unvisited_points = num_voxel - num_iso_endpoint;
t_num_cc = num_iso_endpoint;
while num_unvisited_points > 0
    % Find the starting voxel list index in the voxel list
    for t_start_point_idx = t_start_search_idx : t_num_link_start_voxel
        t_current_id = l_link_start_voxel_idx(t_start_point_idx);
        if voxel_unvisited(t_current_id)
            t_start_search_idx = t_start_point_idx + 1;
            t_num_cc_voxel = 0;   
            break;
        end
    end
    keep_tracking = true;
    while keep_tracking
        keep_tracking = false;
        % Add the current voxel to the connected component voxel list 
        t_num_cc_voxel = t_num_cc_voxel + 1;
        % MATLAB can extend the length of the list automatically if
        % t_num_cc_voxel is larger than the initialized length. 
        t_link_idx(t_num_cc_voxel) = t_current_id;
        voxel_unvisited(t_current_id) = false;
        % Get the neighbors of the current voxel and pick the ONE hasn't
        % been visited. 
        t_neighbor_idx = neighbor_voxel_idx(:, t_current_id);
        for tmp_idx = 1 : 2
            tmp_id = t_neighbor_idx(tmp_idx);
            if tmp_id > 0 % Only one neighbor is unvisited. If the first item is 0, then break the loop directly
                if voxel_unvisited(tmp_id)
                    t_next_id = tmp_id;
                    keep_tracking = true;
                    t_current_id = t_next_id;
                    break;
                end
            else 
                break;
            end
        end
    end
    t_num_cc = t_num_cc + 1;
    num_unvisited_points = num_unvisited_points - t_num_cc_voxel;
    cc.PixelIdxList{t_num_cc} = link_idx_list(t_link_idx(1:t_num_cc_voxel));
end
cc.PixelIdxList = cc.PixelIdxList(1:t_num_cc);
cc.NumObjects = t_num_cc;
cc.Connectivity = 26;
cc.ImageSize = mask_size;
end
