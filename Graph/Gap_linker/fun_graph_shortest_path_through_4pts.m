function connect_link_str = fun_graph_shortest_path_through_4pts(metric, ep1_sub, bp1_sub, bp2_sub, ep2_sub, connectivity)
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

if nargin < 6
    connectivity = 26;
end

local_bbox_size = size(metric);
ep1_ind = sub2ind(local_bbox_size, ep1_sub(1), ep1_sub(2), ep1_sub(3));
bp1_ind = sub2ind(local_bbox_size, bp1_sub(1), bp1_sub(2), bp1_sub(3));
ep2_ind = sub2ind(local_bbox_size, ep2_sub(1), ep2_sub(2), ep2_sub(3));
bp2_ind = sub2ind(local_bbox_size, bp2_sub(1), bp2_sub(2), bp2_sub(3));

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

linker_ep2bv = shortestpath(local_graph, ep1_ind, bp1_ind, 'Method', 'positive');
tmp_s = [];
tmp_t = [];
if numel(linker_ep2bv) > 1
    tmp_s = linker_ep2bv(1:end-1);
    tmp_t = linker_ep2bv(2:end);
end
linker_bv2ep = shortestpath(local_graph, bp2_ind, ep2_ind, 'Method', 'positive');
if numel(linker_bv2ep) > 1
    tmp_s = cat(2, tmp_s, linker_bv2ep(1:end-1));
    tmp_t = cat(2, tmp_t, linker_bv2ep(2:end));
end
if ~isempty(tmp_s)
    local_graph = rmedge(local_graph, tmp_s, tmp_t);
end

linker_bv2bv = shortestpath(local_graph, bp1_ind, bp2_ind, 'Method', 'positive');

connect_link_str.link_ind_w_ep = cat(2, linker_ep2bv, linker_bv2bv(2:end-1), linker_bv2ep)';
connect_link_str.link_sub_w_ep = fun_ind2sub(local_bbox_size, connect_link_str.link_ind_w_ep);
connect_link_str.num_target = 1;
connect_link_str.connectivity = connectivity;
connect_link_str.total_dist = sum(metric(connect_link_str.link_ind_w_ep));
end
