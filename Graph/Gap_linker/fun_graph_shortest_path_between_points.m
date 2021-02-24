function connect_link_str = fun_graph_shortest_path_between_points(metric, point_sub_1, point_sub_2, connectivity)
% fun_graph_shortest_path_between_points finds the shortest path between
% two point, with distance measured by the metric array
% Input: 
%   metric: 3D numerical array, specify the distance between neighboring
%   points in the space
%   point_sub_1: 3-by-1 numerical array, subscript of point 1
%   point_sub_2: N-by-3 numerical array, subscript of the target
%   neighbor_type: 8, 18 or 26. Number of neighbors for searching the
%   shortest path
% Output: 
%   connect_link_str: structure with fields: 
%       link_ind_w_ep: cells of numerical vector, indices of link voxels in the metric array,
%       including the two endpoints
%       total_dist: numerical vector, sum of the metric value of the found link 
%       voxel
%       link_sub_wo_ep: cells of N-by-3 numerical array, subscript of the link
%       voxels, without the subscript of the two input endpoints
%       num_target: number of target
% Author: Xiang Ji
% Date: 01, 09, 2019

if nargin < 4
    connectivity = 26;
end

local_bbox_size = size(metric);
ep_1_local_ind = sub2ind(local_bbox_size, point_sub_1(1), point_sub_1(2), point_sub_1(3));
ep_2_local_ind = sub2ind(local_bbox_size, point_sub_2(:, 1), point_sub_2(:, 2), point_sub_2(:, 3));

% Neighbor list of all the valid indices in the local image
local_valid_mask = ~isinf(metric);
local_valid_ind = find(local_valid_mask); 
metric_pad = padarray(metric, [1, 1, 1], 'both');
local_valid_mask_pad = padarray(local_valid_mask, [1,1,1], 'both');
local_neighbor_ind_add_pad = fun_skeleton_neighbor_add_coeff_3D(local_bbox_size + 2, connectivity, true)';
local_neighbor_ind_add = fun_skeleton_neighbor_add_coeff_3D(local_bbox_size , connectivity, true)';

tmp_sub = fun_ind2sub(local_bbox_size, local_valid_ind);
local_valid_ind_pad = sub2ind(local_bbox_size+2, tmp_sub(:,1) + 1, tmp_sub(:,2) + 1, tmp_sub(:,3) + 1);

local_num_valid_ind = numel(local_valid_ind_pad);
% Construct edge pairs and weight for the edge
local_neighbor_ind_c = cell(local_num_valid_ind, 1);
local_neighbor_weight_c = cell(local_num_valid_ind, 1);
for iter_voxel = 1 : local_num_valid_ind
    tmp_neighbor_ind_pad = local_valid_ind_pad(iter_voxel) +  local_neighbor_ind_add_pad;
    tmp_neighbor_ind = local_valid_ind(iter_voxel) +  local_neighbor_ind_add;
    tmp_neighbor_ind_valid = local_valid_mask_pad(tmp_neighbor_ind_pad);
    local_neighbor_ind_c{iter_voxel} = tmp_neighbor_ind(tmp_neighbor_ind_valid)';    
    local_neighbor_weight_c{iter_voxel} = (metric_pad(local_valid_ind_pad(iter_voxel)) + metric_pad(tmp_neighbor_ind_pad(tmp_neighbor_ind_valid))')./2;
end
local_num_valid_neighbor = cellfun(@numel, local_neighbor_ind_c);

graph_s = repelem(local_valid_ind, local_num_valid_neighbor);
graph_t = cat(1, local_neighbor_ind_c{:});
graph_weight = cat(1, local_neighbor_weight_c{:});
% Search for the shortest path 
assert(all(graph_weight>=0), 'Exist edge of negative weight');
local_graph = graph(graph_s, graph_t, graph_weight);

num_target = size(point_sub_2,1);
if num_target > 1
    connect_link_str.link_ind_w_ep = cell(num_target, 1);
    connect_link_str.link_sub_w_ep = cell(num_target, 1);
    connect_link_str.total_dist = zeros(num_target, 1);
    for iter_target = 1 : num_target
        [connect_link_str.link_ind_w_ep{iter_target} , connect_link_str.total_dist(iter_target)] = ...
            shortestpath(local_graph, ep_1_local_ind, ep_2_local_ind(iter_target), 'Method', 'positive');
%         connect_link_str.link_ind_w_ep{iter_target} = local_valid_ind(path_node_idx);
        connect_link_str.link_sub_w_ep{iter_target} = fun_ind2sub(local_bbox_size, connect_link_str.link_ind_w_ep{iter_target});
    end
    connect_link_str.connectivity = connectivity;
    connect_link_str.num_target = num_target;
else
    [connect_link_str.link_ind_w_ep, connect_link_str.total_dist] = ...
        shortestpath(local_graph, ep_1_local_ind, ep_2_local_ind, 'Method', 'positive');
%     connect_link_str.link_ind_w_ep = local_valid_ind(path_node_idx);
    connect_link_str.link_sub_w_ep = fun_ind2sub(local_bbox_size, connect_link_str.link_ind_w_ep);
    connect_link_str.connectivity = connectivity;
    connect_link_str.num_target = num_target;
end    
end
