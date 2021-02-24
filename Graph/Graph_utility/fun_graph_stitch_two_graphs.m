function graph_c = fun_graph_stitch_two_graphs(graph_1, graph_2, check_overlapQ)
% This function stitch two graph output by fun_skeleton_to_graph to a
% larger graph by integrating information from both graph
% Input: 
%   graph_1, graph_2: two grpah structure defined in fun_skeleton_to_graph
%   check_overlapQ: compute the range of the overlapping part of the link
%   skeleton in the overlapping region. If the range is smaller than
%   min_overlap_range, raise a warning.
% Output: 
%   graph_c: the combined graph
%
% Written by Xiang Ji on September 1st, 2018
% To be updated: 
% 1. Order-preserving link connected component stitching
% 2. Use sparse matrix for storing radius information


if nargin < 3
    check_overlapQ = true;
    min_overlap_range = 20;
end

%% Compute some useful data
graph_1.node = fun_graph_compute_cc_voxel_list_idx(graph_1.node);
graph_1.link = fun_graph_compute_cc_voxel_list_idx(graph_1.link);
graph_2.node = fun_graph_compute_cc_voxel_list_idx(graph_2.node);
graph_2.link = fun_graph_compute_cc_voxel_list_idx(graph_2.link);
%% Determine the stitching direction automatically.
% Assuming only need to stitch in one direction.
bbox_1_global = graph_1.info.bbox_mmxx;
bbox_2_global = graph_2.info.bbox_mmxx;
stitch_direction = unique([find(bbox_1_global(1:3) ~= bbox_2_global(1:3)), ...
    find(bbox_1_global(4:6) ~= bbox_2_global(4:6))]);
assert(isscalar(stitch_direction), 'Need to stitch in more than 1 direction. Have not been implemented yet');
% Make sure that graph_1 always has smaller row/column/layer indices.
if bbox_1_global(stitch_direction) > bbox_2_global(stitch_direction)
    disp('Exchange graph 1 and graph 2.');
    tmp = graph_1;
    graph_1 = graph_2;
    graph_2 = tmp;
    clear tmp;
end
bbox_1_global = graph_1.info.bbox_mmxx;
bbox_2_global = graph_2.info.bbox_mmxx;

% Create a coordinate for these two blocks
pos_common = struct;
% Bounding box of the combined blcok in the global coordinate
pos_common.bbox_global_mmxx = [min(bbox_1_global(1:3), bbox_2_global(1:3)),...
    max(bbox_1_global(4:6), bbox_2_global(4:6))];

pos_common.bbox_global_mmll = [pos_common.bbox_global_mmxx(1:3),...
    pos_common.bbox_global_mmxx(4:6) - pos_common.bbox_global_mmxx(1:3) + 1];

pos_common.bbox_size = pos_common.bbox_global_mmll(4:6);

% Bounding box for block 1 in the common coordinate
pos_common.bbox_block_1_mmxx = [bbox_1_global(1:3) - pos_common.bbox_global_mmxx(1:3),...
    bbox_1_global(4:6) - pos_common.bbox_global_mmxx(1:3)] + 1;

pos_common.bbox_block_1_mmll = [bbox_1_global(1:3) - pos_common.bbox_global_mmxx(1:3),...
    bbox_1_global(4:6) - bbox_1_global(1:3)] + 1;

% Bounding box for block 2 in the common coordinate
pos_common.bbox_block_2_mmxx = [bbox_2_global(1:3) - pos_common.bbox_global_mmxx(1:3),...
    bbox_2_global(4:6) - pos_common.bbox_global_mmxx(1:3)] + 1;

pos_common.bbox_block_2_mmll = [bbox_2_global(1:3) - pos_common.bbox_global_mmxx(1:3),...
    bbox_2_global(4:6) - bbox_2_global(1:3)] + 1;

pos_common.bbox_OR_mmxx  =[pos_common.bbox_block_2_mmxx(1:3), pos_common.bbox_block_1_mmxx(4:6)];
% Divide the common bounding box at the middle in the stitch direction.
pos_common.bbox_mmll = [1,1,1,pos_common.bbox_size];
pos_common.bbox_half_1_mmxx = pos_common.bbox_mmll;
pos_common.bbox_half_1_mmxx(stitch_direction + 3) = round(mean(pos_common.bbox_OR_mmxx([stitch_direction, stitch_direction+3])));


pos_common.num.blk_vox = prod(pos_common.bbox_size);
pos_common.num.bbox_size_pad = pos_common.bbox_size + 2;
pos_common.num.blk_vox_pad = prod(pos_common.num.bbox_size_pad);
%% Convert the linear voxel coordinate in field of pos_ind
% Isopoints
[graph_1.isopoint.c_pos_ind, graph_1.isopoint.c_pos_sub] = fun_linear_index_coordiante_transform3D(graph_1.isopoint.pos_ind, pos_common.bbox_block_1_mmll, [1,1,1,pos_common.bbox_size]);
[graph_2.isopoint.c_pos_ind, graph_2.isopoint.c_pos_sub] = fun_linear_index_coordiante_transform3D(graph_2.isopoint.pos_ind, pos_common.bbox_block_2_mmll, [1,1,1,pos_common.bbox_size]);
% endpoints
[graph_1.endpoint.c_pos_ind, graph_1.endpoint.c_pos_sub] = fun_linear_index_coordiante_transform3D(graph_1.endpoint.pos_ind, pos_common.bbox_block_1_mmll, [1,1,1,pos_common.bbox_size]);
[graph_2.endpoint.c_pos_ind, graph_2.endpoint.c_pos_sub] = fun_linear_index_coordiante_transform3D(graph_2.endpoint.pos_ind, pos_common.bbox_block_2_mmll, [1,1,1,pos_common.bbox_size]);
% Nodes
[graph_1.node.c_pos_ind, graph_1.node.c_pos_sub] = fun_linear_index_coordiante_transform3D(graph_1.node.pos_ind, pos_common.bbox_block_1_mmll, [1,1,1,pos_common.bbox_size]);
[graph_2.node.c_pos_ind, graph_2.node.c_pos_sub] = fun_linear_index_coordiante_transform3D(graph_2.node.pos_ind, pos_common.bbox_block_2_mmll, [1,1,1,pos_common.bbox_size]);
% Links
[graph_1.link.c_pos_ind, graph_1.link.c_pos_sub] = fun_linear_index_coordiante_transform3D(graph_1.link.pos_ind, pos_common.bbox_block_1_mmll, [1,1,1,pos_common.bbox_size]);
[graph_2.link.c_pos_ind, graph_2.link.c_pos_sub] = fun_linear_index_coordiante_transform3D(graph_2.link.pos_ind, pos_common.bbox_block_2_mmll, [1,1,1,pos_common.bbox_size]);

% Spare look up table: common coordinate index to node labels
graph_1.node.map_c_ind_2_label = sparse(graph_1.node.c_pos_ind, ones(graph_1.node.num_voxel,1), ...
    graph_1.node.label, pos_common.num.blk_vox,1);
graph_2.node.map_c_ind_2_label = sparse(graph_2.node.c_pos_ind, ones(graph_2.node.num_voxel,1), ...
    graph_2.node.label, pos_common.num.blk_vox,1);

% Map voxel index in combined mask to voxel index in the original mask 
graph_1.link.map_c_ind_2_ind = sparse(graph_1.link.c_pos_ind, ones(graph_1.link.num_voxel,1),...
    graph_1.link.pos_ind, pos_common.num.blk_vox,1);
graph_2.link.map_c_ind_2_ind = sparse(graph_2.link.c_pos_ind, ones(graph_2.link.num_voxel,1),...
    graph_2.link.pos_ind, pos_common.num.blk_vox,1);

% Get voxel index for each link cc in graph
graph_1.link.cc_c_ind = cell(graph_1.link.num_cc,1);
for tmp_idx = 1 : graph_1.link.num_cc
    graph_1.link.cc_c_ind{tmp_idx} = graph_1.link.c_pos_ind(graph_1.link.cc_idx_range(tmp_idx,1):...
        graph_1.link.cc_idx_range(tmp_idx,2));
end

graph_2.link.cc_c_ind = cell(graph_2.link.num_cc,1);
for tmp_idx = 1 : graph_2.link.num_cc
    graph_2.link.cc_c_ind{tmp_idx} = graph_2.link.c_pos_ind(graph_2.link.cc_idx_range(tmp_idx,1):...
        graph_2.link.cc_idx_range(tmp_idx,2));
end
% Sparse look up table: combined coordinate index to link labels. 
graph_1.link.map_c_ind_2_label = sparse(graph_1.link.c_pos_ind, ones(graph_1.link.num_voxel,1), ...
    graph_1.link.label, pos_common.num.blk_vox,1);
graph_2.link.map_c_ind_2_label = sparse(graph_2.link.c_pos_ind, ones(graph_2.link.num_voxel,1), ...
    graph_2.link.label, pos_common.num.blk_vox,1);
%% Check the range of the overlapping region.
if check_overlapQ
    graph_1.link.c_pos_in_OR_Q = all(bsxfun(@ge, graph_1.link.c_pos_sub, pos_common.bbox_OR_mmxx(1:3)) &...
        bsxfun(@le, graph_1.link.c_pos_sub, pos_common.bbox_OR_mmxx(4:6)),2);
    graph_2.link.c_pos_in_OR_Q = all(bsxfun(@ge, graph_2.link.c_pos_sub, pos_common.bbox_OR_mmxx(1:3)) &...
        bsxfun(@le, graph_2.link.c_pos_sub, pos_common.bbox_OR_mmxx(4:6)),2);
    % Find the link voxels that are not overlapping in two graph, but in the
    % overlapping region
    tmp_c_pos_in_NP_Q = (~full(graph_2.link.map_c_ind_2_label(graph_1.link.c_pos_ind))) & (graph_1.link.c_pos_in_OR_Q);
    tmp_g1_pos_in_NP_sub1 = graph_1.link.c_pos_sub(tmp_c_pos_in_NP_Q, stitch_direction);

    tmp_c_pos_in_NP_Q = (~full(graph_1.link.map_c_ind_2_label(graph_2.link.c_pos_ind))) & (graph_2.link.c_pos_in_OR_Q);
    tmp_g2_pos_in_NP_sub1 = graph_2.link.c_pos_sub(tmp_c_pos_in_NP_Q, stitch_direction);

    tmp_pos1_sort = sort([tmp_g2_pos_in_NP_sub1;tmp_g1_pos_in_NP_sub1], 'ascend');
    [tmp_max_diff, tmp_max_idx] = max(diff(tmp_pos1_sort));
    if tmp_max_diff < min_overlap_range
        warning('The spacing between the cluster of nonoverlapping links is less than 50')
    end
    pos_common.bbox_half_1_mmxx(stitch_direction + 3) = round((tmp_pos1_sort(tmp_max_idx) + tmp_pos1_sort(tmp_max_idx+1))/2);    
end
pos_common.bbox_half_2_mmxx = pos_common.bbox_mmll;
pos_common.bbox_half_2_mmxx(stitch_direction) = pos_common.bbox_half_1_mmxx(stitch_direction + 3) + 1;
%% Find the node for stitching
% Assuming there exist an overlapping region in the middle of the
% overlapping region
% 1. Nodes in graph 1 that have voxels in half 1
% 2. Nodes in graph 2 that have voxels completely in half 2
graph_1.node.c_in_H1_Q = all(bsxfun(@ge, graph_1.node.c_pos_sub, pos_common.bbox_half_1_mmxx(1:3)) &...
    bsxfun(@le, graph_1.node.c_pos_sub, pos_common.bbox_half_1_mmxx(4:6)),2);
graph_2.node.c_in_H1_Q = all(bsxfun(@ge, graph_2.node.c_pos_sub, pos_common.bbox_half_1_mmxx(1:3)) &...
    bsxfun(@le, graph_2.node.c_pos_sub, pos_common.bbox_half_1_mmxx(4:6)),2);

graph_1.node.label_in_H1 = unique(full(graph_1.node.map_c_ind_2_label(graph_1.node.c_pos_ind(graph_1.node.c_in_H1_Q))));
graph_1.node.label_in_H2 = setdiff(graph_1.node.label, graph_1.node.label_in_H1);
graph_2.node.label_in_H1 = unique(full(graph_2.node.map_c_ind_2_label(graph_2.node.c_pos_ind(graph_2.node.c_in_H1_Q))));
graph_2.node.label_in_H2 = setdiff(graph_2.node.label, graph_2.node.label_in_H1);
%% Stitch node
graph_1.node.cc_c_ind = cell(graph_1.node.num_cc,1);
for tmp_idx = 1 : graph_1.node.num_cc
    graph_1.node.cc_c_ind{tmp_idx} = graph_1.node.c_pos_ind(graph_1.node.cc_idx_range(tmp_idx,1):...
        graph_1.node.cc_idx_range(tmp_idx,2));
end
graph_2.node.cc_c_ind = cell(graph_2.node.num_cc,1);
for tmp_idx = 1 : graph_2.node.num_cc
    graph_2.node.cc_c_ind{tmp_idx} = graph_2.node.c_pos_ind(graph_2.node.cc_idx_range(tmp_idx,1):...
        graph_2.node.cc_idx_range(tmp_idx,2));
end
graph_c = struct;
graph_c.num.blk_vox = pos_common.num.blk_vox;
graph_c.node.idx_last_from_1 = length(graph_1.node.cc_c_ind(graph_1.node.label_in_H1));
graph_c.node.cc_ind = cat(1, graph_1.node.cc_c_ind(graph_1.node.label_in_H1),...
    graph_2.node.cc_c_ind(graph_2.node.label_in_H2));

if isfield(graph_1.node, 'radius') && isfield(graph_2.node, 'radius')
    graph_c.node.radius = cat(1, graph_1.node.radius(cat(2, graph_1.node.cc_list_idx{graph_1.node.label_in_H1})),...
        graph_2.node.radius(cat(2, graph_2.node.cc_list_idx{graph_2.node.label_in_H2})));
end

graph_c.node.num_cc = length(graph_c.node.cc_ind);
graph_c.node.pos_ind = cat(1, graph_c.node.cc_ind{:});
graph_c.node.num_voxel = length(graph_c.node.pos_ind);
graph_c.node.num_voxel_per_cc = cellfun(@length, graph_c.node.cc_ind);
graph_c.node.label = repelem(1:graph_c.node.num_cc, graph_c.node.num_voxel_per_cc)';
graph_c.node.map_ind_2_label = sparse(graph_c.node.pos_ind, ...
    ones(graph_c.node.num_voxel,1), ...
    graph_c.node.label, graph_c.num.blk_vox,1);
graph_c.node.connected_link_label = cell(graph_c.node.num_cc,1);
graph_c.node.num_link = zeros(graph_c.node.num_cc,1);
%% Stitch Links
%% Find overlapping link voxels
% Find the link voxels in the overlapping region
graph_1.link.c_in_H1_Q = all(bsxfun(@ge, graph_1.link.c_pos_sub, pos_common.bbox_half_1_mmxx(1:3)) &...
    bsxfun(@le, graph_1.link.c_pos_sub, pos_common.bbox_half_1_mmxx(4:6)),2);
graph_2.link.c_in_H1_Q = all(bsxfun(@ge, graph_2.link.c_pos_sub, pos_common.bbox_half_1_mmxx(1:3)) &...
    bsxfun(@le, graph_2.link.c_pos_sub, pos_common.bbox_half_1_mmxx(4:6)),2);
%% Classify links
% Graph 1:
tmp_link_label_not_all_in_H1 = unique(graph_1.link.label(~graph_1.link.c_in_H1_Q));
% Label for links not completely out of B1S ( May be completely in B1S )
tmp_label_in_H1 = unique(graph_1.link.label(graph_1.link.c_in_H1_Q));
% Label for links completely in B1S
graph_1.link.label_all_in_H1 = setdiff(tmp_label_in_H1, tmp_link_label_not_all_in_H1);
% Label for links cross the boundary of B1S
graph_1.link.label_part_in_H1 = intersect(tmp_label_in_H1, tmp_link_label_not_all_in_H1);

% graph 2
tmp_link_label_not_all_in_H1 = unique(graph_2.link.label(~graph_2.link.c_in_H1_Q));
% Label for links not completely out of H1 ( May be completely in B1S )
tmp_label_in_H1 = unique(graph_2.link.label(graph_2.link.c_in_H1_Q));
% Label for links cross the boundary of H2
graph_2.link.label_part_in_H2  = intersect(tmp_label_in_H1, tmp_link_label_not_all_in_H1);
% Label for links completely in H2
graph_2.link.label_all_in_H2 = setdiff(tmp_link_label_not_all_in_H1, graph_2.link.label_part_in_H2);
assert(length(graph_1.link.label_part_in_H1) == length(graph_2.link.label_part_in_H2),...
    'Number of link crossing the middle plane of the overlapping region is not the same in two graph.');
%% graph_c.link initialization
graph_c.link = struct;
est_num_link = length(graph_1.link.label_all_in_H1) + length(graph_2.link.label_all_in_H2) + ...
    length(graph_2.link.label_part_in_H2) + 1000;
graph_c.link.cc_ind = cell(est_num_link,1);
graph_c.link.num_node = zeros(est_num_link,1);
graph_c.link.connected_node_label = zeros(est_num_link,2);
graph_c.link.num_cc = 0;
if isfield(graph_1.link, 'radius') && isfield(graph_2.link, 'radius')
    need_to_stitch_radius_Q = true;
    graph_c.link.radius = cell(3,1);
end
%% Stitch all the links that are completely in the first half
% disp('Stitch links that are completely in the first half');
tmp_num_label = length(graph_1.link.label_all_in_H1);
graph_c.link.cc_ind(graph_c.link.num_cc + 1:graph_c.link.num_cc + tmp_num_label ) = graph_1.link.cc_c_ind(graph_1.link.label_all_in_H1);
tmp_num_node_list = graph_1.link.num_node(graph_1.link.label_all_in_H1);
tmp_connected_node_label_list = graph_1.link.connected_node_label(graph_1.link.label_all_in_H1,:);
for tmp_idx = 1 : tmp_num_label
    tmp_link_label_c = graph_c.link.num_cc + tmp_idx;
    tmp_connected_node_label_in_1 = tmp_connected_node_label_list(tmp_idx,1:tmp_num_node_list(tmp_idx));
    
    if ~isempty(tmp_connected_node_label_in_1) % If this link connects to any nodes
        for tmp_node_label = tmp_connected_node_label_in_1
            % Position in graph 1
            tmp_node_c_ind = graph_1.node.cc_c_ind{tmp_node_label};
            % Node label in graph c
            tmp_node_label_in_c = unique(full(graph_c.node.map_ind_2_label(tmp_node_c_ind)));
            % Need to deal with link connected to nodes in graph 1 but not
            % graph 2
            if nnz(tmp_node_label_in_c)>0
                tmp_node_label_in_c = tmp_node_label_in_c(tmp_node_label_in_c>0);
                assert(length(tmp_node_label_in_c)<2, 'Node in graph 1 correspond to multiple nodes in the combined graph');
                % Record information to node
                graph_c.node.num_link(tmp_node_label_in_c) = graph_c.node.num_link(tmp_node_label_in_c) + 1;
                graph_c.node.connected_link_label{tmp_node_label_in_c}(graph_c.node.num_link(tmp_node_label_in_c)) =  tmp_link_label_c;
                % Add information to link
                graph_c.link.num_node(tmp_link_label_c) = graph_c.link.num_node(tmp_link_label_c) + 1;
                graph_c.link.connected_node_label(tmp_link_label_c, graph_c.link.num_node(tmp_link_label_c)) = tmp_node_label_in_c;
            else
                warning('Node in the combined graph is deleted. Link label %d\n', tmp_link_label_c);
                % Link connected to node that has been deleted.
                % Find the endpint and add to endpoint list
            end
        end
    end
end
if need_to_stitch_radius_Q
    graph_c.link.radius{1} = graph_1.link.radius(cat(2, graph_1.link.cc_list_idx{graph_1.link.label_all_in_H1}));
end
graph_c.link.num_cc = graph_c.link.num_cc + tmp_num_label;
%% Modify the voxel list for links corss the boundary and link nodes
% Get B1S part of the voxel indices list from graph_1
% disp('Stitch links that cross the boundary of the first and the second half.');
tmp_num_label = length(graph_2.link.label_part_in_H2);
tmp_radius_cell_array = cell(tmp_num_label, 1);
for tmp_idx = 1 : tmp_num_label
    tmp_link_label_c = graph_c.link.num_cc + tmp_idx;
    % Link label in graph 2
    tmp_label_2 = graph_2.link.label_part_in_H2(tmp_idx);
    % Link voxel list idx in grpah 2
    tmp_list_idx = graph_2.link.cc_list_idx{tmp_label_2};
    % Link voxel ind in common coordinate
    tmp_ind = graph_2.link.cc_c_ind{tmp_label_2};
    % Link label in graph 1 that connected with
    tmp_label_in_BS1 = full(graph_1.link.map_c_ind_2_label(tmp_ind));
    tmp_label_in_BS1 = unique(tmp_label_in_BS1(tmp_label_in_BS1>0));
    % Connected link list idx in graph 1
    tmp_idx_in_B1 = cat(2,graph_1.link.cc_list_idx{tmp_label_in_BS1});
    % Part of the voxels in graph 1 that are in H1
    tmp_pos_1_in_B1S_Q = all( bsxfun(@le, graph_1.link.c_pos_sub(tmp_idx_in_B1,:),...
        pos_common.bbox_half_1_mmxx(4:6)),2);
    tmp_cc_ind_B1S = graph_1.link.c_pos_ind(tmp_idx_in_B1(tmp_pos_1_in_B1S_Q));
    % Part of the voxels in graph 2 that are in B2S
    tmp_pos_2_in_B2S_Q = all( bsxfun(@ge, graph_2.link.c_pos_sub(tmp_list_idx,:),...
        pos_common.bbox_half_2_mmxx(1:3)),2);
    tmp_cc_ind_B2S = graph_2.link.c_pos_ind(tmp_list_idx(tmp_pos_2_in_B2S_Q));
    % Combined link voxel list and add to graph_c
    graph_c.link.cc_ind{tmp_link_label_c} = cat(1,tmp_cc_ind_B1S, tmp_cc_ind_B2S);
    if need_to_stitch_radius_Q
        tmp_radius_B1S = graph_1.link.radius(tmp_idx_in_B1(tmp_pos_1_in_B1S_Q));
        tmp_radius_B2S = graph_2.link.radius(tmp_list_idx(tmp_pos_2_in_B2S_Q));
        tmp_radius_cell_array{tmp_idx} = cat(1,tmp_radius_B1S, tmp_radius_B2S);
    end
    % Node connected to in graph 1
    tmp_node_label = graph_1.link.connected_node_label(tmp_label_in_BS1,:);
    tmp_node_label = tmp_node_label(tmp_node_label>0);
    % Check if the connected nodes are in B1S. If yes, connect them.
    %     tmp_node_c_ind = zeros(length(tmp_node_label(:)),1);
    for tmp_node_idx = 1:length(tmp_node_label(:))
        tmp_node_list_idx = graph_1.node.cc_list_idx{tmp_node_label(tmp_node_idx)};
        tmp_node_in_B1SQ = all( bsxfun(@le, graph_1.node.c_pos_sub(tmp_node_list_idx,:),...
            pos_common.bbox_half_1_mmxx(4:6)),2);
        tmp_node_ind = graph_1.node.c_pos_ind(tmp_node_list_idx(tmp_node_in_B1SQ));
        if ~isempty(tmp_node_ind)
            tmp_node_label_in_c = unique(full(graph_c.node.map_ind_2_label(tmp_node_ind)));
            assert(numel(tmp_node_label_in_c)<2, 'Node in graph 1 correspond to multiple node in the combined grpah');
            %             fprintf('Add node %d from graph 1\n', tmp_node_label_in_c);
            % Add link to node
            graph_c.node.num_link(tmp_node_label_in_c) = graph_c.node.num_link(tmp_node_label_in_c) + 1;
            graph_c.node.connected_link_label{tmp_node_label_in_c}(graph_c.node.num_link(tmp_node_label_in_c)) =  tmp_link_label_c;
            % Add information to link
            graph_c.link.num_node(tmp_link_label_c) = graph_c.link.num_node(tmp_link_label_c) + 1;
            graph_c.link.connected_node_label(tmp_link_label_c, graph_c.link.num_node(tmp_link_label_c)) = tmp_node_label_in_c;
        end
    end
    
    % Node connected to in graph 2
    tmp_node_label = graph_2.link.connected_node_label(tmp_label_2,:);
    tmp_node_label = tmp_node_label(tmp_node_label>0);
    % Check if the connected nodes are in B1S. If yes, connect them.
    %     tmp_node_c_ind = zeros(length(tmp_node_label(:)));
    for tmp_node_idx = 1:length(tmp_node_label)
        tmp_node_list_idx = graph_2.node.cc_list_idx{tmp_node_label(tmp_node_idx)};
        tmp_node_in_B2SQ = all( bsxfun(@ge, graph_2.node.c_pos_sub(tmp_node_list_idx,:),...
            pos_common.bbox_half_2_mmxx(1:3)),2);
        tmp_node_ind = graph_2.node.c_pos_ind(tmp_node_list_idx(tmp_node_in_B2SQ));
        if ~isempty(tmp_node_ind)
            tmp_node_label_in_c = unique(full(graph_c.node.map_ind_2_label(tmp_node_ind)));
            assert(numel(tmp_node_label_in_c)<2, 'Node in graph 1 correspond to multiple node in the combined grpah');
            %             fprintf('Add node %d from graph 2\n', tmp_node_label_in_c);
            % Add link to node
            graph_c.node.num_link(tmp_node_label_in_c) = graph_c.node.num_link(tmp_node_label_in_c) + 1;
            graph_c.node.connected_link_label{tmp_node_label_in_c}(graph_c.node.num_link(tmp_node_label_in_c)) =  tmp_link_label_c;
            % Add information to link
            graph_c.link.num_node(tmp_link_label_c) = graph_c.link.num_node(tmp_link_label_c) + 1;
            graph_c.link.connected_node_label(tmp_link_label_c, graph_c.link.num_node(tmp_link_label_c)) = tmp_node_label_in_c;
        end
    end
    %     check_link_voxel_1 = setdiff(cat(1,graph_c.link.cc_ind{:}), combined_graph_gt.link.pos_ind);
    %     assert(isempty(check_link_voxel_1),'Introduce wrong voxels');
end
if need_to_stitch_radius_Q
    graph_c.link.radius{2} = cat(1, tmp_radius_cell_array{:});
end
graph_c.link.num_cc = graph_c.link.num_cc + tmp_num_label;
%% Stitch all the links that are completely in B2S
% disp('Stitch links that are completely in B2S');
tmp_num_label = length(graph_2.link.label_all_in_H2);
graph_c.link.cc_ind(graph_c.link.num_cc+ 1:graph_c.link.num_cc + tmp_num_label ) = graph_2.link.cc_c_ind(graph_2.link.label_all_in_H2);
tmp_num_node = graph_2.link.num_node(graph_2.link.label_all_in_H2);
tmp_connected_node_label = graph_2.link.connected_node_label(graph_2.link.label_all_in_H2,:);
for tmp_idx = 1 : tmp_num_label
    tmp_link_label_c = graph_c.link.num_cc + tmp_idx;
    tmp_connected_node_label_in_1 = tmp_connected_node_label(tmp_idx,1:tmp_num_node(tmp_idx));
    if ~isempty(tmp_connected_node_label_in_1)
        for tmp_node_label = tmp_connected_node_label_in_1
            tmp_node_c_ind = graph_2.node.cc_c_ind{tmp_node_label};
            tmp_node_label_in_c = unique(full(graph_c.node.map_ind_2_label(tmp_node_c_ind)));
            if nnz(tmp_node_label_in_c)>0
                tmp_node_label_in_c = tmp_node_label_in_c(tmp_node_label_in_c>0);
                assert(length(tmp_node_label_in_c)<2, 'Node in graph 1 correspond to multiple nodes in the combined graph');
                % Record information to node
                graph_c.node.num_link(tmp_node_label_in_c) = graph_c.node.num_link(tmp_node_label_in_c) + 1;
                graph_c.node.connected_link_label{tmp_node_label_in_c}(graph_c.node.num_link(tmp_node_label_in_c)) =  tmp_link_label_c;
                % Add information to link
                graph_c.link.num_node(tmp_link_label_c) = graph_c.link.num_node(tmp_link_label_c) + 1;
                graph_c.link.connected_node_label(tmp_link_label_c, graph_c.link.num_node(tmp_link_label_c)) = tmp_node_label_in_c;
            else
                warning('Node in the combined graph is deleted. Link label %d\n', tmp_link_label_c);
                % Link connected to node that has been deleted.
                % Find the endpint and add to endpoint list
            end
        end
    end
end
graph_c.link.num_cc = graph_c.link.num_cc + tmp_num_label;

graph_c.link.pos_ind = cat(1, graph_c.link.cc_ind{1:graph_c.link.num_cc});
graph_c.link.num_voxel = length(graph_c.link.pos_ind);
graph_c.link.num_voxel_per_cc = cellfun(@length, graph_c.link.cc_ind(1:graph_c.link.num_cc));
graph_c.link.label = repelem(1:graph_c.link.num_cc, graph_c.link.num_voxel_per_cc(1:graph_c.link.num_cc));
graph_c.link.map_ind_2_label = sparse(graph_c.link.pos_ind, ones(graph_c.link.num_voxel,1), ...
    graph_c.link.label, pos_common.num.blk_vox,1);
if need_to_stitch_radius_Q
    graph_c.link.radius{3} = graph_2.link.radius(cat(2, graph_2.link.cc_list_idx{graph_2.link.label_all_in_H2}));
    graph_c.link.radius = cat(1, graph_c.link.radius{:});
end
%% Endpoints and isolated points.
graph_1.endpoint.c_in_H1_Q = all(bsxfun(@ge, graph_1.endpoint.c_pos_sub, pos_common.bbox_half_1_mmxx(1:3)) &...
    bsxfun(@le, graph_1.endpoint.c_pos_sub, pos_common.bbox_half_1_mmxx(4:6)),2);
graph_2.endpoint.c_in_H2_Q = all(bsxfun(@ge, graph_2.endpoint.c_pos_sub, pos_common.bbox_half_2_mmxx(1:3)) &...
    bsxfun(@le, graph_2.endpoint.c_pos_sub, pos_common.bbox_half_2_mmxx(4:6)),2);

graph_c.endpoint.num_voxel = nnz(graph_1.endpoint.c_in_H1_Q) + nnz(graph_2.endpoint.c_in_H2_Q);
graph_c.endpoint.pos_ind = cat(1, graph_1.endpoint.c_pos_ind(graph_1.endpoint.c_in_H1_Q), ...
    graph_2.endpoint.c_pos_ind(graph_2.endpoint.c_in_H2_Q));
if need_to_stitch_radius_Q
    graph_c.endpoint.radius = cat(1, graph_1.endpoint.radius(graph_1.endpoint.c_in_H1_Q), ...
    graph_2.endpoint.radius(graph_2.endpoint.c_in_H2_Q));
end
graph_c.endpoint.link_label = full(graph_c.link.map_ind_2_label(graph_c.endpoint.pos_ind));
assert(~ismember(0, graph_c.endpoint.link_label), 'Exist endpoints do not belong to any link');

graph_1.isopoint.c_in_H1_Q = all(bsxfun(@ge, graph_1.isopoint.c_pos_sub, pos_common.bbox_half_1_mmxx(1:3)) &...
    bsxfun(@le, graph_1.isopoint.c_pos_sub, pos_common.bbox_half_1_mmxx(4:6)),2);
graph_2.isopoint.c_in_H2_Q = all(bsxfun(@ge, graph_2.isopoint.c_pos_sub, pos_common.bbox_half_2_mmxx(1:3)) &...
    bsxfun(@le, graph_2.isopoint.c_pos_sub, pos_common.bbox_half_2_mmxx(4:6)),2);
graph_c.isopoint.num_voxel = nnz(graph_1.isopoint.c_in_H1_Q) + nnz(graph_2.isopoint.c_in_H2_Q);
graph_c.isopoint.pos_ind = cat(1, graph_1.isopoint.c_pos_ind(graph_1.isopoint.c_in_H1_Q), ...
    graph_2.isopoint.c_pos_ind(graph_2.isopoint.c_in_H2_Q));
if need_to_stitch_radius_Q
    graph_c.isopoint.radius = cat(1, graph_1.isopoint.radius(graph_1.isopoint.c_in_H1_Q), ...
    graph_2.isopoint.radius(graph_2.isopoint.c_in_H2_Q));
end
%% Graph information and parameters
graph_c.info.dataset_name = graph_1.info.dataset_name;
graph_c.info.stack = graph_1.info.stack;
graph_c.info.layer = union(graph_1.info.layer, graph_2.info.layer);
graph_c.info.idx_1 = union(graph_1.info.idx_1, graph_2.info.idx_1);
graph_c.info.idx_2 = union(graph_1.info.idx_2, graph_2.info.idx_2);
graph_c.info.bbox_mmxx = pos_common.bbox_global_mmxx;
graph_c.info.bbox_mmll = pos_common.bbox_global_mmll;

graph_c.num.mask_size = graph_c.info.bbox_mmll(4:6);
graph_c.num.skeleton_voxel = graph_c.link.num_voxel + graph_c.node.num_voxel + ...
    graph_c.endpoint.num_voxel + graph_c.isopoint.num_voxel;
graph_c.num.block_voxel = prod(graph_c.num.mask_size);
%% Verification
% %% Generate 'ground truth'
% % Stitch direction dependent.
% combined_cell_size = [1,1,1];
% combined_cell_size(stitch_direction) = 2;
% combined_block_mask = cell(combined_cell_size);
% combined_block_mask{1} = segment_info_1.distance_transform > 0;
% combined_block_mask{2} = segment_info_2.distance_transform > 0;
% [combined_block_mask, tmp]= fun_stitch_blocks_with_overlap( combined_block_mask, 240, 1);
% combined_skeleton_gt = bwskel(combined_block_mask);
% combined_graph_gt = fun_skeleton_to_graph_v2(find(combined_skeleton_gt), size(combined_skeleton_gt));
% %% Check with the ground truth
% tmp_error_node_gt_but_c = setdiff(combined_graph_gt.node.pos_ind, graph_c.node.pos_ind);
% tmp_error_node_c_but_gt = setdiff(graph_c.node.pos_ind, combined_graph_gt.node.pos_ind);
% if isempty(tmp_error_node_c_but_gt)
%     disp('All the stitched nodes are the nodes in ground truth');
% else
%     disp('Exist stitched nodes not existing in ground truth');
% end
% if isempty(tmp_error_node_gt_but_c)
%     disp('All the ground truth nodes are in the stitched graph');
% else
%     disp('Exist ground truth nodes not existing the stithced graph');
%     disp(tmp_error_node_gt_but_c(:)');
% end
% %% Check link and connectivity with ground truth
% check_link_voxel_1 = setdiff(cat(1,graph_c.link.cc_ind{:}), combined_graph_gt.link.pos_ind);
% if numel(check_link_voxel_1) == 0
%     fprintf('All the link voxels in stitched graph are in ground truth\n');
% else
%     fprintf('%d link voxels in stitched graph are not in ground truth\n', length(check_link_voxel_1));
% end
% % check_link_voxel_2 = setdiff(combined_graph_gt.link.pos_ind, link_voxel_before_stitch);
% check_link_voxel_2 = setdiff(combined_graph_gt.link.pos_ind, cat(1,graph_c.link.cc_ind{:}));
% if numel(check_link_voxel_2) == 0
%     fprintf('All the link voxels in ground truth are in the stitched graph\n');
% else
%     fprintf('%d link voxels in ground truth are not in the stitched graph\n', numel(check_link_voxel_2));
% end
%
% if length(combined_graph_gt.link.pos_ind) ~= length(cat(1,graph_c.link.cc_ind{:}))
%     fprintf('The number of the link voxels in the ground truth is the same as the stitched graph.\n');
% end
% %% For each node, check all the link voxel connected to it
% disp('Check if all the node connected to the correct set of voxel');
% for tmp_idx = 1 : graph_c.node.num_cc
%     % Node position in graph_c
%     node_ind = graph_c.node.cc_ind{tmp_idx};
%     % Node label in gt
%     gt_node_label = unique(full(combined_graph_gt.node.map_ind_2_label(node_ind)));
%     assert(isscalar(gt_node_label), 'gt_node_label is not a scalar');
%     % link voxel of this node in grpah_c
%     link_voxel_list_c = cat(1, graph_c.link.cc_ind{graph_c.node.connected_link_label{tmp_idx}});
%     % link voxel of this node in gt
%     link_voxel_list_gt = cat(1, combined_graph_gt.link.cc_ind{combined_graph_gt.node.connected_link_label{gt_node_label}});
%     if ~isempty(setdiff(link_voxel_list_c,link_voxel_list_gt))
%         error('Node %d is incorrect. It has link voxels not in ground truth.\n', tmp_idx);
%     end
%     if ~isempty(setdiff(link_voxel_list_gt,link_voxel_list_c))
%         error('Node %d is incorrect. Ground truth voxels are missing\n', tmp_idx);
%     end
% end
% disp('Check done. All the links connected to the same nodes.');
% %% For each link, check the location of all the node it attached to
% disp('Check if all the link connected to the correct node');
% for tmp_idx = 1 : graph_c.link.num_cc
%     node_label = graph_c.link.connected_node_label(tmp_idx,1:graph_c.link.num_node(tmp_idx));
%     node_ind = cat(1, graph_c.node.cc_ind{node_label});
%     node_label_gt = unique(full(combined_graph_gt.node.map_ind_2_label(node_ind)));
%     if ismember(0, node_label_gt)
%         error('Incorrect node label %d. Exist node voxel not in ground truth\n', tmp_idx);
%     end
%     node_label_gt = node_label_gt(node_label_gt>0);
%     node_ind_gt = cat(1, combined_graph_gt.node.cc_ind{node_label_gt});
%     if ~isempty(setdiff(node_ind_gt,node_ind))
%         error('Link %d is incorrect. It has node voxels not in ground truth.\n', tmp_idx);
%     end
%     if ~isempty(setdiff(node_ind, node_ind_gt))
%         error('Link %d is incorrect. Ground truth voxels are missing\n', tmp_idx);
%     end
% end
% disp('Check done. All the node connected to the same links.');
end

function input_str = fun_graph_compute_cc_voxel_list_idx(input_str)
% This function compute extra useful information for later processing. 
% Input:
%   input_str: graph structure is generated by fun_skeleton_to_graph_v2
% Output:
%   input_str: add three more fields to the input structure

input_str.cc_idx_range = zeros(input_str.num_cc,2);
% input_str.cc_c_ind = cell(input_str.num_cc,1);
input_str.cc_idx_range(:,2) = cumsum(input_str.num_voxel_per_cc(:));
input_str.cc_idx_range(:,1) = input_str.cc_idx_range(:,2) - ...
    input_str.num_voxel_per_cc(:) + 1;

input_str.cc_list_idx = cell(input_str.num_cc,1);
for tmp_idx = 1 : input_str.num_cc
    input_str.cc_list_idx{tmp_idx} = input_str.cc_idx_range(tmp_idx,1) : ...
        input_str.cc_idx_range(tmp_idx,2);
%     input_str.cc_c_ind{tmp_idx} = input_str.c_pos_ind(input_str.cc_idx_range(tmp_idx,1):...
%         input_str.cc_idx_range(tmp_idx,2));
end
end

function [local_linear_ind_t , vargout]= fun_linear_index_coordiante_transform3D(local_linear_ind, bbox_cur_in_target_mmll, bbox_target_mmll)
% This function convert the array indices in one block to the array indices
% in the other block. 
% Input: 
%     local_linear_ind: N-by-1 double precision array

%     bbox_cur_in_target_mmll: the bounding box from which
%     local_linear_ind is subtracted, [min_pos1, min_pos2, min_pos3, l1,
%     l2, l3]. min_pos1 is the minimum index in the first dimensiton (row)
%     in the coordinate of the target

%     bbox_target_mmll: normally equals [1,1,1, size_of_bbox]

if any(bbox_cur_in_target_mmll(1:3) < bbox_target_mmll(1:3) | ...
        (bbox_cur_in_target_mmll(1:3) + bbox_cur_in_target_mmll(4:6)) > ...
        (bbox_target_mmll(1:3) + bbox_target_mmll(4:6)))
    error('The current bouding box is out of the range of the target bounding box');
end

[pos1, pos2, pos3] = ind2sub(bbox_cur_in_target_mmll(4:6), local_linear_ind);
pos1_t = pos1 + bbox_cur_in_target_mmll(1) - bbox_target_mmll(1);
pos2_t = pos2 + bbox_cur_in_target_mmll(2) - bbox_target_mmll(2);
pos3_t = pos3 + bbox_cur_in_target_mmll(3) - bbox_target_mmll(3);
local_linear_ind_t = sub2ind(bbox_target_mmll(4:6), pos1_t, pos2_t, pos3_t);

if nargout > 1
    vargout = cat(2, pos1_t, pos2_t, pos3_t);
end

end