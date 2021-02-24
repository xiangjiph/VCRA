function [input_graph, new_link_info] = fun_graph_pruning_convert_node_to_link(input_graph, node_label)
% fun_graph_pruning_convert_node_to_link is part of the graph pruning
% algorithm. It merges the node voxels with the connected link (preserve
% order) and update graph structure (without relabeling).
% Input: 
%   input_graph: graph structure created by fun_skeleton_to_graph
%   node_label: scalar, the label of the node to be converted. 
% Output: 
%   input_graph: updated graph structure
%   new_link_info: structure with fileds:
%       label: new link label
%       length: length of the new link 
%       num_endpoint: number of endpoints of this link
% Written by Xiang Ji on Nov 11, 2018
% Modified by Xiang Ji on Jan 27, 2019
% 1. Fix a bug in determine neighboring ind
% Modified by Xiang Ji on APr 26, 2019
% 1. Fix a bug caused by the updated fun_skeleton_to_graph, where the
% number of edges each node connect to is counted by link label, so
% self-loops are only counted once. 

% Find the remaining two links
if ~isfield(input_graph, 'tmp_link_ind_old_to_new')
    input_graph.tmp_link_ind_old_to_new = 0 : input_graph.link.num_cc;
end
if ~isscalar(node_label)
    error('This function can only handle one node each time');
end
link_label = input_graph.node.connected_link_label{node_label};
link_label = input_graph.tmp_link_ind_old_to_new(link_label + 1);
num_link_label = numel(link_label);
if  num_link_label > 2
    error('More than 2 links connected to the node');
elseif num_link_label < 2
    error('Debug');
end
link_cc_ind_1 = input_graph.link.cc_ind{link_label(1)};
link_cc_ind_2 = input_graph.link.cc_ind{link_label(2)};
link_length_1 = numel(link_cc_ind_1);
link_length_2 = numel(link_cc_ind_2);
if link_label(1) == link_label(2)
%     fprintf('Node %d is part of an isolated loop\n', node_label);
    new_link_info.is_loop = true;
    new_link_info.label = link_label(1);
    new_link_info.num_endpoint = 0;
    new_link_info.length = link_length_1;
    return;
end
% After merging two links, two node need to be modified - The two
% node can be endpoint...
connected_node_label_1 = input_graph.link.connected_node_label(link_label(1),:);
connected_node_label_1 = connected_node_label_1(connected_node_label_1 ~= node_label);
connected_node_label_2 = input_graph.link.connected_node_label(link_label(2),:);
connected_node_label_2 = connected_node_label_2(connected_node_label_2 ~= node_label);


if ~isempty(connected_node_label_1) && ~isempty(connected_node_label_2)
    % do nothing
elseif isempty(connected_node_label_2) && ~isempty(connected_node_label_1)
    new_link_info.is_loop = true;
    new_link_info.label = link_label(1);
    new_link_info.num_endpoint = 0;
    new_link_info.length = link_length_1;
    return
elseif isempty(connected_node_label_1) && ~isempty(connected_node_label_2)
    new_link_info.is_loop = true;
    new_link_info.label = link_label(2);
    new_link_info.num_endpoint = 0;
    new_link_info.length = link_length_2;
    return
else
    new_link_info.is_loop = true;
    new_link_info.label = nan;
    new_link_info.num_endpoint = 0;
    new_link_info.length = nan;
%     fprintf('The node probably connected to two self-links. Debug...\n');
    return    
end
    

% For debug
% fprintf('Number of endpoint in the link: %d\n', numel(connected_node_label == 0));
% Neigoboring link voxel of the node:
node_ind = input_graph.node.cc_ind{node_label};
node_size = numel(node_ind);
% Find the neighboring link labels for each node voxel
tmp_sub = fun_ind2sub(input_graph.num.mask_size, node_ind);
switch size(tmp_sub, 2)
    case 3
        node_ind_pad = sub2ind(input_graph.num.mask_size_pad, tmp_sub(:,1)+1, ...
            tmp_sub(:,2)+1, tmp_sub(:,3)+1);
    case 2
        node_ind_pad = sub2ind(input_graph.num.mask_size_pad, tmp_sub(:,1)+1, ...
            tmp_sub(:,2)+1);
end
node_n26_ind_pad = bsxfun(@plus, node_ind_pad', input_graph.num.neighbor_add_pad);
tmp_sub_pad = fun_ind2sub(input_graph.num.mask_size_pad, node_n26_ind_pad(:));
num_dim = size(tmp_sub_pad, 2);
% Restrict the subscript between 2 - mask size + 1
tmp_sub_valid_Q = all( bsxfun(@ge, tmp_sub_pad, 2 .* ones(1, num_dim)), 2) & ...
    all( bsxfun(@le, tmp_sub_pad, input_graph.num.mask_size_pad - 1), 2);

node_n26_ind = bsxfun(@plus, node_ind', input_graph.num.neighbor_add);
node_n26_link_label = zeros(size(node_n26_ind));
node_n26_link_label(tmp_sub_valid_Q) = full(input_graph.link.map_ind_2_label(node_n26_ind(tmp_sub_valid_Q)));
% Map the old link label to a new one ( the map_ind_to_label is not updated
% due to performance cost) 
node_n26_link_label = input_graph.tmp_link_ind_old_to_new(node_n26_link_label + 1);
if node_size == 1
%     node_n26_ind = input_graph.node.cc_ind{node_label} + input_graph.num.neighbor_add;
%     % Rescrict the indices within the bounding box
%     node_n26_ind = node_n26_ind(node_n26_ind>0 & node_n26_ind<= input_graph.num.block_voxel);
%     node_n26_link_label = full(input_graph.link.map_ind_2_label(node_n26_ind)); 
%     % Conver the old label to the new label
%     node_n26_link_label = input_graph.tmp_link_ind_old_to_new(node_n26_link_label + 1);
%     node_n26_ind_validQ = node_n26_link_label > 0;
%     node_n26_ind = node_n26_ind(node_n26_ind_validQ);
%     node_n26_link_label = node_n26_link_label(node_n26_ind_validQ);
    % Determine the order of the cc ind. Connect two link to the
    % old node voxel.
    connected_ind_1 = node_n26_ind(node_n26_link_label == link_label(1));
    connected_ind_2 = node_n26_ind(node_n26_link_label == link_label(2));
    if link_cc_ind_1(1) == connected_ind_1
        link_cc_ind_1 = flip(link_cc_ind_1);
    elseif link_cc_ind_1(link_length_1) ~= connected_ind_1
        error('The node is not the endpoint of the link cc.');
    end
    if link_cc_ind_2(link_length_2) == connected_ind_2
        link_cc_ind_2 = flip(link_cc_ind_2);
    elseif link_cc_ind_2(1) ~= connected_ind_2
        error('The node is not the endpoint of the link cc.');
    end
    node_voxel_in_link = node_ind;
else
    % 1. Find the link voxels in two link cc that need to be
    % connected.
    % 2. Use the two link voxels, as well as all of the node voxels
    % to constuct the connectivity matrix. Use the relative
    % position to define the edge values
    % 3. Find the shortest path starting from one link voxel to
    % another one. Modify the graph accordingly.
    
    % It is impossible to have a node voxel connected to two voxels
    % of the same link cc.
    [tmp_neighbor_idx, tmp_voxel_idx, tmp_neighbor_link_label] = find(node_n26_link_label);
    tmp_neighbor_ind = zeros(size(tmp_neighbor_idx));
    for tmp_idx = 1 : numel(tmp_neighbor_idx)
        tmp_neighbor_ind(tmp_idx) = node_n26_ind(tmp_neighbor_idx(tmp_idx), tmp_voxel_idx(tmp_idx));
    end
    
    connected_idx_1 = find(tmp_neighbor_link_label == link_label(1));
    connected_ind_1 = tmp_neighbor_ind(connected_idx_1);
    connected_idx_1 = tmp_voxel_idx(connected_idx_1);
    if numel(connected_ind_1) > 1
%         disp('More than 1 node voxels connect to the same link voxel');
        connected_ind_1 = connected_ind_1(1);
        connected_idx_1 = connected_idx_1(1);
    end
    connected_idx_2 = find(tmp_neighbor_link_label == link_label(2));
    connected_ind_2 = tmp_neighbor_ind(connected_idx_2);
    connected_idx_2 = tmp_voxel_idx(connected_idx_2);
    if numel(connected_ind_2) > 1
%         disp('More than 1 node voxels connect to the same link voxel');
        connected_ind_2 = connected_ind_2(1);
        connected_idx_2 = connected_idx_2(1);
    end
    
    if link_cc_ind_1(1) == connected_ind_1
        link_cc_ind_1 = flip(link_cc_ind_1);
    elseif link_cc_ind_1(link_length_1) ~= connected_ind_1
        error('The node is not the endpoint of the link cc.');
    end
    if link_cc_ind_2(link_length_2) == connected_ind_2
        link_cc_ind_2 = flip(link_cc_ind_2);
    elseif link_cc_ind_2(1) ~= connected_ind_2
        error('The node is not the endpoint of the link cc.');
    end
    
    % Construct connectivity matrix - Treat 26 neighbors the same
    tmp_cc_ind = [connected_ind_1;node_ind;connected_ind_2];
    tmp_cc_numel = numel(tmp_cc_ind);
    tmp_connectivity_matrix = zeros(tmp_cc_numel);
    tmp_connectivity_matrix(1+connected_idx_1, 1) = 1 + connected_idx_1;
    tmp_connectivity_matrix(1, 1+connected_idx_1) = 1;
    tmp_connectivity_matrix(1+connected_idx_2, tmp_cc_numel) = 1+connected_idx_2;
    tmp_connectivity_matrix(tmp_cc_numel, 1+connected_idx_2) = tmp_cc_numel;
    % Internal node connectivity
    for tmp_idx_1 = 1 : node_size
        for tmp_idx_2 = tmp_idx_1+1 : node_size
            if (tmp_idx_1 ~= tmp_idx_2) && any(node_n26_ind(:, tmp_idx_1) == node_ind(tmp_idx_2))
                tmp_connectivity_matrix(tmp_idx_1+1, tmp_idx_2+1) = tmp_idx_1 + 1;
                tmp_connectivity_matrix(tmp_idx_2+1, tmp_idx_1+1) = tmp_idx_2 + 1;
            end
        end
    end
    tmp_connectivity_matrix = tmp_connectivity_matrix + tmp_connectivity_matrix';
    % Evoking the shortest path finding function is time consuming. Deal
    % with the simple case first    
    if tmp_connectivity_matrix(1 + connected_idx_1, tmp_cc_numel)
        % Single voxel connected to both links
        node_voxel_in_link = node_ind(connected_idx_1);
    elseif tmp_connectivity_matrix(1 + connected_idx_1, 1 + connected_idx_2)
        node_voxel_in_link = node_ind([connected_idx_1; connected_idx_2]);
    elseif node_size == 3
        if connected_idx_1 == 1
            if connected_idx_2 == 2
                middle_idx = 3;
            else
                middle_idx = 2;
            end
        elseif connected_idx_1 == 2
            if connected_idx_2 == 1
                middle_idx = 3;
            else
                middle_idx = 1;
            end
        else
            if connected_idx_2 == 1
                middle_idx = 2;
            else
                middle_idx = 1;
            end
        end
        node_voxel_in_link = node_ind([connected_idx_1; middle_idx; connected_idx_2]);
    else
%         warning('Node of size larger than 3');
        tmp_cc_graph = graph(tmp_connectivity_matrix>0);
        tmp_path = shortestpath(tmp_cc_graph, 1, tmp_cc_numel, 'Method', 'unweighted');
        node_voxel_in_link = tmp_cc_ind(tmp_path(2:end-1));
    end
%     % For debug
%     tmp_cc_graph = graph(tmp_connectivity_matrix>0);
%     tmp_path = shortestpath(tmp_cc_graph, 1, tmp_cc_numel, 'Method', 'unweighted');
%     node_voxel_in_link_1 = tmp_cc_ind(tmp_path(2:end-1));
%     if any(node_voxel_in_link_1 ~= node_voxel_in_link)
%         error('The simplified version of shortestpath gives the wrong answer');
%     end
end
new_link_cc = [link_cc_ind_1; node_voxel_in_link; link_cc_ind_2];

% Update Node information: modify the connected link label for the
% second node
if connected_node_label_1 ~= 0 || (connected_node_label_1 == 0 && connected_node_label_2 == 0)
    % If link 1 connect to a node
    kept_node_label = connected_node_label_1;
    modified_node_label = connected_node_label_2;
    keep_link_label = link_label(1);
    modified_link_label = link_label(2);
%     modified_link_ind = [node_voxel_in_link; link_cc_ind_2];
else 
    % if link 1 connect to an endpoint and link 2 connect to a node
    kept_node_label = connected_node_label_2;
%     modified_link_ind = [link_cc_ind_1; node_voxel_in_link];
    modified_node_label = 0;
    keep_link_label = link_label(2);
    modified_link_label = link_label(1);
end

if (modified_node_label == 0)
    % if two endpoints, change the link label of endpoint 2 to the link
    % label of endpoint 1;
    % if only one endpoints, change the link label of the endpoint to the
    % new link label.
    input_graph.endpoint.link_label(input_graph.endpoint.link_label == modified_link_label) = keep_link_label;
elseif kept_node_label ~= 0 % and modified_node_label ~= 0
    % connect to two nodes
    link_label_node_mod =  input_graph.node.connected_link_label{modified_node_label};
    link_label_node_mod(link_label_node_mod == modified_link_label) = keep_link_label;
    input_graph.node.connected_link_label{modified_node_label} = link_label_node_mod;
end
% Update link information: Update the new connected component
input_graph.link.cc_ind{keep_link_label} = new_link_cc;
% input_graph.link.map_ind_2_label(modified_link_ind) = keep_link_label;
% input_graph.tmp_link_ind_old_to_new(modified_link_label + 1) = keep_link_label;
% The modified link can be previously modified and therefore all the
% original links that have been merged into the current modified_link_label
% should be modified. 
input_graph.tmp_link_ind_old_to_new(input_graph.tmp_link_ind_old_to_new == modified_link_label) = keep_link_label;
% input_graph.link.map_ind_2_label(node_voxel_in_link) = keep_link_label;
% Delete link 2
input_graph.link.cc_ind{modified_link_label} = [];
input_graph.link.connected_node_label(modified_link_label,:) = 0;
% Modify the node connected to
input_graph.link.connected_node_label(keep_link_label,:) = [kept_node_label, modified_node_label];
% Modify the number of voxels in the link
new_link_length = length(new_link_cc);
input_graph.link.num_voxel_per_cc(keep_link_label) = new_link_length;
% Remove node
input_graph.node.cc_ind{node_label} = [];
% input_graph.node.map_ind_2_label(node_ind) = 0;

new_link_info.label = keep_link_label;
new_link_info.length = new_link_length;
new_link_info.is_loop = false;
% Possible to create loop......
if modified_node_label ~= 0
    new_link_info.num_endpoint = 0;
elseif kept_node_label ~= 0
    new_link_info.num_endpoint = 1;
else
    new_link_info.num_endpoint = 2;
end
% Deal with loop...
if kept_node_label == modified_node_label
    new_link_info.is_loop = true;
end
    
end