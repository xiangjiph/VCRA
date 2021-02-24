function graph_str = fun_skeleton_to_graph(voxel_list, mask_size)
% This function computes the graph from the skeleton. 
% Input:
%   voxel_list: N-by-1 double precision array of skeleton voxel indices
%   mask_size: 3-by-1 double precision array. size of the mask(skeleton
%   array)
%   min_link_length: links with the number of voxel(including the endpoint)
%   smaller than min_link_length are not recorded.. 
% Output: graph_str with fields:
%   node:structure with fields:
%         pos_ind: linear index position of the node voxel in the mask 
%         num_voxel, num_cc(#node): each node can have more than 1 voxels
%         label: N-by-1 double precision array. specify the label of the
%         node voxel in pos_ind;
%         cc_ind: cell array, each contains the linear index of the voxel
%         in each node
%         num_link: number of link that this node attached
%         connected_link_label: labels of links that this node joins.
%         map_ind_2_label: N-by-1 sparse double precision matrix for
%         mapping the voxel linear index to its label.
%   num: structure with fields
%         mask_size(3-by-1), skeleton_voxel(number of skeleton voxel),
%         mask_size_pad( = mask_size +2, i.e. pad one layer of 0 on both
%         size), blk_vol(number of voxels in the block, or mask),
%         neighbor_add_pad(26-by-1 double precision array for finding the
%         26 neighbors, in the padded array), and block_voxel_pad.
%   link, isopoint and endpoint: strucutres with similar fileds with node. 
%
% Written by Xiang Ji on Aug 23, 2018
% 1. Faster than original VIDA implementation, with better strucutred output
% 2. View each node cluster as single node. 
% 3. Output voxel list in link are not in order(in terms of vessel direction) 
% 4. Deal with skeleton voxels on the boundary 
%
% Modified by Xiang Ji on Aug 28, 2018
% 1. Change the order of pos_ind, which is now extracted from
% cc.PixelIdxList. 
% 2. The cell array storing list index of the voxels in the pos_ind can be
% easily calculated since pos_ind is the concatenate of PixelIdxList. For
% saving storage space, this cell array is not added to the structures. 
% Plan:
% 1. If the mex implementation of finding connected componennt on sparse
% matrix is compariable to the MATLAB implementation, the step for
% calculating the list voxel index of each connected component can be
% removed by modifying the information recorded in fun_cc_in_sparse_matrix.
%
% Modified by Xiang Ji on Sep 10, 2018
% 1. Vectorize the process for finding the number of neighbors
% 2. Add minimum link length threshod. 
%
% Modified by Xiang Ji on Sep 27, 2018
% 1. Use the neighborhood label of each voxel to track links (replace
% fun_cc_in_sparse). The current implementation is faster than
% reconstructing skeleton and use bwconcomp. Moreover, the collected
% voxels starts from a endpoint/node and end at the endpoint/node. 
% 2. Delete the link lengh thresholding. The previous implementation might
% delete short linking voxels between two nodes. Here, we
% don't want to delete any voxels from the skeleton. 
%
% Modified by Xiang Ji on Oct 29, 2018
% 1. Fix a bug on finding link connected component.
% 2. Deal with isolated loops
% Modified by Xiang Ji on Nov 7, 2018
% 1. Allow input 3D skeleton logical array
% 
% Modified by Xiang Ji on Apr 25, 2019
% 1. Links that form loops and connect to the same node have two identical,
% nonzero connected_to_node_label value. Previously, one of the label is 0.
% The value in the node structure ramain the same. Original implementation
% is saved as fun_skeleton_to_graph_v3
% 
% Modified by Xiang Ji on Aug 13, 2019
% 1. Track along the isolated loop connected component, instead of just
% finding the connected components. 
% The original implementation is saved as fun_skeleton_to_graph_v4
%% Initialization
if ~isvector(voxel_list)
    % If the input is the skeleton (3D logical array), convert it to voxel
    % list. 
    if nargin < 2
        mask_size = size(voxel_list);
    end
    voxel_list = find(voxel_list);
elseif ~isa(voxel_list, 'double')
    voxel_list = double(voxel_list);    
end
num.mask_size = mask_size;
num.mask_size_pad = num.mask_size + 2;
num.block_voxel = prod(num.mask_size);
num.block_voxel_pad = prod(num.mask_size_pad);
% 26 neighbors relative indices position array:
tmp1 = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3];
tmp2 = [1 1 1 2 2 2 3 3 3 1 1 1 2 2 2 3 3 3 1 1 1 2 2 2 3 3 3];
tmp3 = [1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3]; 
num.neighbor_add_pad = sub2ind(num.mask_size_pad, tmp1(:), tmp2(:), tmp3(:));
num.neighbor_add_pad = num.neighbor_add_pad - num.neighbor_add_pad(14);
num.neighbor_add_pad(14) = [];

num.neighbor_add = sub2ind(num.mask_size, tmp1(:), tmp2(:), tmp3(:));
num.neighbor_add = num.neighbor_add - num.neighbor_add(14);
num.neighbor_add(14) = [];
%% Generate 1D sparse matrix representation of the skeleton array
[pos_1, pos_2, pos_3] = ind2sub(num.mask_size, voxel_list);
num.skeleton_voxel = numel(voxel_list);
voxel_idx_padded = sub2ind(num.mask_size_pad, pos_1 + 1, pos_2 + 1, pos_3 + 1); 
sp_skl = sparse(voxel_idx_padded, ones(num.skeleton_voxel, 1), ...
    1:num.skeleton_voxel, num.block_voxel_pad, 1);
clear pos_1 pos_2 pos_3
%% Classify voxels
% list of neighbor voxel index in voxel_list for each voxel. 
l_voxel_neighbor_idx = full(sp_skl(bsxfun(@plus, voxel_idx_padded', num.neighbor_add_pad)));
% Sort the neighbor for link tracking in the next section
l_voxel_neighbor_idx = sort(l_voxel_neighbor_idx, 1, 'descend');
l_num_neighbor = sum(l_voxel_neighbor_idx>0,1);
% List of endpoint index in the voxel list
l_ep_idx = find(l_num_neighbor == 1);
% List of node index in the voxel list
l_nd_idx = find(l_num_neighbor > 2);

l_isop_Q = (l_num_neighbor == 0);
%% Track along links
l_neighbor_label_of_nd_idx = l_voxel_neighbor_idx(:, l_nd_idx);
l_neighbor_label_of_nd_idx = l_neighbor_label_of_nd_idx(l_neighbor_label_of_nd_idx>0);
l_neighbor_label_of_nd_idx = unique(l_neighbor_label_of_nd_idx);
% The following vector contains both the endpoint indices and link point
% indices. 
l_neighbor_link_label_of_nd_idx = setdiff(l_neighbor_label_of_nd_idx, l_nd_idx);
% The start tacking points are the union of the link points in the neighbor
% to the node points, or the end points ( in case that the links have two
% endpoints)
l_link_start_voxel_idx = union(l_neighbor_link_label_of_nd_idx, l_ep_idx);
% Label the voxels need to visit
voxel_unvisited = true(num.skeleton_voxel,1);
voxel_unvisited(l_nd_idx) = false;
voxel_unvisited(l_isop_Q) = false;

[link_cc, num_unvisited_points, voxel_unvisited] = fun_skeleton_get_link_cc(...
    voxel_list, l_voxel_neighbor_idx, l_link_start_voxel_idx, voxel_unvisited);
%% Construct graph
graph_str = struct;
graph_str.num = num;

% Link information
link_length_list = cellfun(@length, link_cc.PixelIdxList');
% Should be very carefully when trying to delete some links. Links,
% endpoints and nodes should be modified in the same time. 
graph_str.link.num_cc = link_cc.NumObjects;
graph_str.link.cc_ind = link_cc.PixelIdxList';
graph_str.link.pos_ind = cat(1, graph_str.link.cc_ind{:});
graph_str.link.num_voxel = numel(graph_str.link.pos_ind);
graph_str.link.num_voxel_per_cc = link_length_list;

% Transpose to make a column vector
graph_str.link.label = repelem(1:graph_str.link.num_cc, graph_str.link.num_voxel_per_cc)';
graph_str.link.map_ind_2_label = sparse(graph_str.link.pos_ind, ...
    ones(graph_str.link.num_voxel,1), ...
    graph_str.link.label, ...
    graph_str.num.block_voxel,1);

% Endpoint information (Notice that endpoints are a subset of link points)
graph_str.endpoint.pos_ind = voxel_list(l_ep_idx);
graph_str.endpoint.link_label = full(graph_str.link.map_ind_2_label(graph_str.endpoint.pos_ind));
% The following three line of code is used when some links have been
% deleted. 
% endpoint_in_link_Q = graph_str.endpoint.link_label>0;
% graph_str.endpoint.pos_ind = graph_str.endpoint.pos_ind(endpoint_in_link_Q);
% graph_str.endpoint.link_label = graph_str.endpoint.link_label(endpoint_in_link_Q);
graph_str.endpoint.num_voxel = numel(l_ep_idx);
graph_str.endpoint.map_ind_2_label = sparse(graph_str.endpoint.pos_ind, ones(graph_str.endpoint.num_voxel, 1),...
    1 : graph_str.endpoint.num_voxel, graph_str.num.block_voxel, 1);
%% Isolated point
graph_str.isopoint.pos_ind = voxel_list(l_isop_Q);
graph_str.isopoint.num_voxel = nnz(graph_str.isopoint.pos_ind);
%% Isolated loop
if num_unvisited_points > 0
    try
        [loop_cc, num_unvisited_points_1, ~] = fun_skeleton_get_link_cc(voxel_list, l_voxel_neighbor_idx, find(voxel_unvisited), voxel_unvisited);
        assert(num_unvisited_points_1 == 0, 'Still exist unvisited voxels after finding for isolated loop connected component');
    catch ME
        fprintf('The new version of finding isolated loop connected components fail. Use the old version\n');
        fprintf('Error identifier %s\n Error message: %s\n', ME.identifier, ME.message);
        loop_cc = fun_cc_in_sparse_matrix(voxel_list(voxel_unvisited), mask_size);
    end
else
    loop_cc.PixelIdxList = cell(0,1);
end
graph_str.isoloop.cc_ind = loop_cc.PixelIdxList';
graph_str.isoloop.num_cc = numel(graph_str.isoloop.cc_ind);
graph_str.isoloop.pos_ind = cat(1, graph_str.isoloop.cc_ind{:});
%% Correcting node information
graph_str.link.connected_node_label = zeros(graph_str.link.num_cc,2);
graph_str.link.num_node = zeros(graph_str.link.num_cc,1);
if isempty(l_nd_idx)
    graph_str.node.num_voxel = 0;
    graph_str.node.num_cc = 0;
    graph_str.node.cc_ind = {};
    graph_str.node.pos_ind = [];
    graph_str.node.num_voxel_per_cc = [];
    graph_str.node.label = [];
    graph_str.node.map_ind_2_label = [];
    graph_str.node.num_link = [];
    graph_str.node.connected_link_label = {};
    return;
end 
node_cc = fun_cc_in_sparse_matrix(voxel_list(l_nd_idx), mask_size);
% tmp_mask = false(num.mask_size);
% tmp_mask(voxel_list(l_nd_idx)) = true;
% node_cc_1 = bwconncomp(tmp_mask);

% Node information should be processed here because node can become link
% voxel if a short link with one endpoint is removed. 
% Node information
graph_str.node.num_voxel = nnz(l_nd_idx);
graph_str.node.num_cc = node_cc.NumObjects;
% Transpose the cell array
graph_str.node.cc_ind = node_cc.PixelIdxList';
graph_str.node.pos_ind = cat(1, node_cc.PixelIdxList{:});
graph_str.node.num_voxel_per_cc = cellfun(@length, graph_str.node.cc_ind);

% Transpose the vector to make it a column vector
graph_str.node.label = repelem(1:graph_str.node.num_cc, graph_str.node.num_voxel_per_cc)';
graph_str.node.map_ind_2_label = sparse(graph_str.node.pos_ind, ...
    ones(graph_str.node.num_voxel,1), ...
    graph_str.node.label, ...
    graph_str.num.block_voxel,1);
%% Connect nodes with links
% Since the order of the voxel indices has been changed, the linear
% position of each voxel cannot be extracted by the logical array anymore.
% Need to re-calculate or find a list indices map. 
[pos_1, pos_2, pos_3] = ind2sub(num.mask_size, graph_str.link.pos_ind);
link_pos_ind_pad = sub2ind(num.mask_size_pad, pos_1 + 1, pos_2 + 1, pos_3 + 1); 
[pos_1, pos_2, pos_3] = ind2sub(num.mask_size, graph_str.node.pos_ind);
node_pos_ind_pad = sub2ind(num.mask_size_pad, pos_1 + 1, pos_2 + 1, pos_3 + 1); 
% Be careful about the construction of the map look up table!!!
link_map_ind_pad_2_lable = sparse(link_pos_ind_pad, ones(graph_str.link.num_voxel, 1), ...
    graph_str.link.label, num.block_voxel_pad,1);

graph_str.node.num_link = zeros(graph_str.node.num_cc,1);
graph_str.node.connected_link_label = cell(graph_str.node.num_cc,1);
% for each voxel in the node list, look at the 26 neighbors in the padded
% link point look up table. 
for tmp_node_idx = 1 : graph_str.node.num_voxel
    tmp_node_label = graph_str.node.label(tmp_node_idx);

    node_neighbor_link_label = full(link_map_ind_pad_2_lable(node_pos_ind_pad(tmp_node_idx) + num.neighbor_add_pad));
    node_neighbor_link_label = node_neighbor_link_label(node_neighbor_link_label>0);
    % Is seems to be possilbe for finding multiple voxels of the same link
    % in the 26 neighbor. (???)
    % One possibility is self-loop
    num_node_neighbor_link_voxels = numel(node_neighbor_link_label);
%     if num_node_neighbor_link_voxels ~= num_neighbor_seg
%         fprintf('Multiple node voxel connected to the same link\n');
%     end
    % Add information to node
    graph_str.node.connected_link_label{tmp_node_label} = cat(1, graph_str.node.connected_link_label{tmp_node_label}, node_neighbor_link_label);
    % Add information to link
    for tmp_neighbor_idx = 1 : num_node_neighbor_link_voxels
        tmp_l_voxel_neighbor_idx = node_neighbor_link_label(tmp_neighbor_idx);
        % Add information to link        
        graph_str.link.num_node(tmp_l_voxel_neighbor_idx) = graph_str.link.num_node(tmp_l_voxel_neighbor_idx) + 1;
        assert(graph_str.link.num_node(tmp_l_voxel_neighbor_idx) < 3, 'The link connect to more than 2 nodes voxels');
        graph_str.link.connected_node_label(tmp_l_voxel_neighbor_idx, graph_str.link.num_node(tmp_l_voxel_neighbor_idx)) = tmp_node_label;
    end
end
graph_str.node.connected_link_label = cellfun(@unique, graph_str.node.connected_link_label, 'UniformOutput', false);
graph_str.node.num_link = cellfun(@numel, graph_str.node.connected_link_label);
% graph_str.link.connected_node_label = sort(graph_str.link.connected_node_label,2, 'descend');
% Completeness check
% if isempty(setdiff(voxel_list, cat(1, graph_str.node.pos_ind, ...
%         graph_str.link.pos_ind, graph_str.isopoint.pos_ind)))
%     disp('All the voxels in the skeleton have been classified');
% else
%     disp('Bug exist');
% end
% if isempty(setdiff(cat(1, graph_str.node.pos_ind, ...
%         graph_str.link.pos_ind, graph_str.isopoint.pos_ind),voxel_list))
%     disp('All the voxels in the graph comes from the skeleton');
% else
%     disp('Bug exist');
% end

end
