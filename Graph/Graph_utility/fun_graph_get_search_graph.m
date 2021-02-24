function search_graph = fun_graph_get_search_graph(metric, connectivity)
% fun_graph_get_search_graph convert the metric to a graph for shortest
% path searching 
% Input: 
%   metric: 3D numerical array. Inf means the voxel is not connected to
%   anything other voxel. The distance between two neighboring voxels are
%   the average of their voxel value in the metric (to make the distance
%   measure symmetric). 
%   connectivity: definition of neighboring voxel, can be 8, 18 or
%   26(default)
% Output: 
%   search_graph: undirected MATLAB graph structure
% 
% Implemented by Xiang Ji on 03/14/2019
%
if nargin < 2
    connectivity = 26;
end

local_bbox_size = size(metric);
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
% Construct graph
assert(all(graph_weight>=0), 'Exist edge of negative weight');
search_graph = graph(graph_s, graph_t, graph_weight);
end