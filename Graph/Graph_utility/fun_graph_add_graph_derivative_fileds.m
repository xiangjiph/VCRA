function input_graph = fun_graph_add_graph_derivative_fileds(input_graph)
% fun_graph_add_graph_derivative_fileds conputes the fields in the graph
% that can be derived from {cc_ind, connected_(node)/(link)_label}
% Input: 
%   input_graph: structure with fields 
%       link: structure with fields cc_ind and connected_node_label
%       node: structure with fields cc_ind and connected_link_label
%       endpoint: structure with field pos_ind
%       * For the definition of these fields, see fun_skeleton_to_graph
% Output: 
%   input_graph: updated structure
%
% Implemented by Xiang Ji on 02/18/2019

if isfield(input_graph, 'link')
    input_graph.link.num_cc = numel(input_graph.link.cc_ind);
    input_graph.link.pos_ind = cat(1, input_graph.link.cc_ind{:});
    input_graph.link.num_voxel = numel(input_graph.link.pos_ind);
    input_graph.link.num_voxel_per_cc = cellfun(@numel, input_graph.link.cc_ind);
    input_graph.link.num_node = sum(input_graph.link.connected_node_label>0,2);
    input_graph.link.label = repelem(1:input_graph.link.num_cc, input_graph.link.num_voxel_per_cc)';
    input_graph.link.map_ind_2_label = sparse(input_graph.link.pos_ind, ...
        ones(input_graph.link.num_voxel,1), ...
        input_graph.link.label, ...
        input_graph.num.block_voxel,1);
end

if isfield(input_graph, 'node')
    input_graph.node.num_cc = numel(input_graph.node.cc_ind);
    input_graph.node.pos_ind = cat(1, input_graph.node.cc_ind{:});
    input_graph.node.num_voxel = numel(input_graph.node.pos_ind);
    input_graph.node.num_voxel_per_cc = cellfun(@numel, input_graph.node.cc_ind);
    input_graph.node.num_link = cellfun(@numel, input_graph.node.connected_link_label);
    input_graph.node.label = repelem(1:input_graph.node.num_cc, input_graph.node.num_voxel_per_cc)';
    input_graph.node.map_ind_2_label = sparse(input_graph.node.pos_ind, ...
        ones(input_graph.node.num_voxel,1), ...
        input_graph.node.label, ...
        input_graph.num.block_voxel,1);
end

if isfield(input_graph, 'endpoint')
    input_graph.endpoint.link_label =  full(input_graph.link.map_ind_2_label(input_graph.endpoint.pos_ind));
    valid_endpoint_voxel_Q = (input_graph.endpoint.link_label > 0);
    input_graph.endpoint.pos_ind = input_graph.endpoint.pos_ind(valid_endpoint_voxel_Q);
    input_graph.endpoint.link_label = input_graph.endpoint.link_label(valid_endpoint_voxel_Q);
    input_graph.endpoint.num_voxel = numel(input_graph.endpoint.pos_ind);
    input_graph.endpoint.map_ind_2_label = sparse(input_graph.endpoint.pos_ind, ones(input_graph.endpoint.num_voxel, 1),...
        1 : input_graph.endpoint.num_voxel, input_graph.num.block_voxel, 1);
end
% input_graph.num.skeleton_voxel = input_graph.link.num_voxel + input_graph.node.num_voxel + input_graph.isopoint.num_voxel;
input_graph.num.skeleton_voxel = input_graph.link.num_voxel + input_graph.node.num_voxel;

end