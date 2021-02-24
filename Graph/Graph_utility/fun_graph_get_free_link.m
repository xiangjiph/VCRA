function link_info = fun_graph_get_free_link(vessel_graph, internal_offset, sortQ)
% 
%
%
%
% Note: 
% Link with two endpoints can appear in the list of link_info.ep1, because
% one of the endpoint can be out of the internal_offset and therefore the
% link effectively has only one internal endpoint
if nargin < 2
    internal_offset = 0;
    sortQ = false;
elseif nargin < 3
    sortQ = false;    
end
block_size = vessel_graph.num.mask_size;
if internal_offset > 0    
    internal_bbox = [ones([1,3])*internal_offset,  block_size - internal_offset + 1];
    [endpoint_pos_sub(:,1), endpoint_pos_sub(:,2), ...
        endpoint_pos_sub(:,3)] = ind2sub(block_size, vessel_graph.endpoint.pos_ind);
    endpoint_internal_Q = all(bsxfun(@gt, endpoint_pos_sub, ...
        internal_bbox(1:3)),2) & all(bsxfun(@lt, endpoint_pos_sub, ...
        internal_bbox(4:6)),2);
    endpoint_ind = vessel_graph.endpoint.pos_ind(endpoint_internal_Q);
else
    endpoint_ind = vessel_graph.endpoint.pos_ind;
end
link_label_ep = full(vessel_graph.link.map_ind_2_label(endpoint_ind));
[link_label_unique, ~, tmp_idx] = unique(link_label_ep, 'stable');
num_internal_ep = cellfun(@numel, fun_bin_data_to_idx_list(tmp_idx));
link_label_ep2 = link_label_unique(num_internal_ep == 2);

link_info = struct;
link_info.num.mask_size = vessel_graph.num.mask_size;
link_info.ep_sub = fun_ind2sub(link_info.num.mask_size, endpoint_ind);
link_info.ep_ind = endpoint_ind;
link_info.num.ep = numel(endpoint_ind);

link_info.ep2.link_label = link_label_ep2;
link_info.ep2.num_cc = numel(link_info.ep2.link_label);
link_info.ep2.link_cc_ind = vessel_graph.link.cc_ind(link_info.ep2.link_label);
link_info.ep2.num_voxel_per_cc = vessel_graph.link.num_voxel_per_cc(link_info.ep2.link_label);
% Endpoints in the same connected component
assert(all(link_info.ep2.num_voxel_per_cc  >= 2), 'Exist link connected component with 1 voxel');
link_info.ep2.ep_pair_ind = zeros(2, link_info.ep2.num_cc);
for iter1 = 1 : link_info.ep2.num_cc
    link_info.ep2.ep_pair_ind(:, iter1) = link_info.ep2.link_cc_ind{iter1}([1; end]);
end
link_info.ep2.ep_ind = link_info.ep2.ep_pair_ind(:);
link_info.ep2.ep_link_label = full(vessel_graph.link.map_ind_2_label(link_info.ep2.ep_ind));
link_info.ep2.ep_sub = fun_ind2sub(link_info.num.mask_size, link_info.ep2.ep_ind);

assert(~any(any(vessel_graph.link.connected_node_label(link_info.ep2.link_label, :),2)), 'Link with two endpoints connected to a node');
link_info.ep2.link_cc_sub = cell(link_info.ep2.num_cc, 1);
link_info.ep2.length = zeros(link_info.ep2.num_cc, 1);
for iter1 = 1 : link_info.ep2.num_cc
    link_info.ep2.link_cc_sub{iter1} = fun_ind2sub(link_info.num.mask_size, link_info.ep2.link_cc_ind{iter1});
    link_info.ep2.length(iter1) = fun_graph_sub_to_length(link_info.ep2.link_cc_sub{iter1}, 1);
end

[link_info.ep1.ep_ind, tmp_idx] = setdiff(endpoint_ind, link_info.ep2.ep_ind);
link_info.ep1.ep_sub = fun_ind2sub(link_info.num.mask_size, link_info.ep1.ep_ind);

link_info.ep1.link_label = link_label_ep(tmp_idx);
link_info.ep1.num_cc = numel(link_info.ep1.link_label);
link_info.ep1.link_cc_ind = vessel_graph.link.cc_ind(link_info.ep1.link_label);
link_info.ep1.num_voxel_per_cc = vessel_graph.link.num_voxel_per_cc(link_info.ep1.link_label);
link_info.ep1.connected_node_label = vessel_graph.link.connected_node_label(link_info.ep1.link_label, :);

link_info.ep1.link_cc_sub = cell(link_info.ep1.num_cc, 1);
link_info.ep1.length = zeros(link_info.ep1.num_cc, 1);
for iter1 = 1 : link_info.ep1.num_cc
    % Set the endpoint to be the last element of the list
    tmp_ind = link_info.ep1.link_cc_ind{iter1};
    if numel(tmp_ind) > 1
        if tmp_ind(1) == link_info.ep1.ep_ind(iter1)
            tmp_ind = flip(tmp_ind);
            link_info.ep1.link_cc_ind{iter1} = tmp_ind;
        end
    end    
    link_info.ep1.link_cc_sub{iter1} = fun_ind2sub(link_info.num.mask_size, tmp_ind);
    link_info.ep1.length(iter1) = fun_graph_sub_to_length(link_info.ep1.link_cc_sub{iter1}, 1);
end

assert(~any(link_info.ep1.connected_node_label(:,2)), 'Link with one endpoint connect to two nodes');
% Link with one endpoint can cannot to no node since here we only consider
% internal endpoint. 
% assert(all(link_info.ep1.connected_node_label(:,1)), 'Link with one endpoint connect to no nodes');
%% Sort link from long to short
if sortQ
    [link_info.ep1.num_voxel_per_cc, tmp_link_idx] = sort(link_info.ep1.num_voxel_per_cc, 'descend');
    tmp_fn = fieldnames(link_info.ep1);
    for iter_fn = 1 : numel(tmp_fn)
        if ~strcmp(tmp_fn{iter_fn} , 'num_cc') && ~strcmp(tmp_fn{iter_fn} , 'num_voxel_per_cc')
            link_info.ep1.(tmp_fn{iter_fn}) = link_info.ep1.(tmp_fn{iter_fn})(tmp_link_idx, :);
        end
    end
end
end