function [input_graph, new_endpoint_ind] = fun_graph_delete_internal_links(input_graph, link_label_to_remove)
% fun_graph_delete_internal_links remove the links (connected to two nodes)
% and return the new graph
% If th_hair_length = 0, do not remove any of the newly created endpoint
th_hair_length = 0;
%% Remove 
num_link_to_remove = numel(link_label_to_remove);
num_new_endpoint = 0;
new_endpoint_ind = zeros(ceil(num_link_to_remove/5), 1);
if num_link_to_remove == 0
%     disp('No link needs to be removed');
    return;
else
%     fprintf('Need to remove %d links\n', num_link_to_remove);
end
input_graph.tmp_link_ind_old_to_new = 0 : input_graph.link.num_cc;
link_to_be_deleted_length = input_graph.link.num_voxel_per_cc(link_label_to_remove);
for int_ep_idx = 1 : num_link_to_remove
    link_label = link_label_to_remove(int_ep_idx);  
    % Map the link label to the updated one.
    link_label = input_graph.tmp_link_ind_old_to_new(link_label + 1);
    if input_graph.link.num_voxel_per_cc(link_label) ~= link_to_be_deleted_length(int_ep_idx)
%         fprintf('The length of the link %d to be removed has changed. Skip this link\n',...
%             link_label);
        continue;
    end
    while ~isempty(link_label) && (link_label ~= 0)  
        % Remove the link:
        [input_graph, updated_node_info] = fun_graph_pruning_remove_single_link(input_graph, link_label, false);
        if isempty(updated_node_info.label)
%             disp('This link does not connect to any node');
            link_label = [];
            continue;
        end
        [~, tmp_idx] = unique(updated_node_info.label);
        for tmp_node_idx = tmp_idx'
            connected_node_label = updated_node_info.label(tmp_node_idx);
            if updated_node_info.degree(tmp_node_idx) == 2
%                 fprintf('Delete node and join two link connected components\n');
                [input_graph, new_link_info] = fun_graph_pruning_convert_node_to_link(input_graph, connected_node_label);
                if (new_link_info.length > th_hair_length) || (new_link_info.num_endpoint == 0)
                    link_label = [];
                else
                    link_label = new_link_info.label;
                end
            elseif updated_node_info.degree(tmp_node_idx) == 1
                % After deleting a loop, it's possible that the originally
                % connected node only connect to one links. One option here
                % is to convert the node into an endpoint and update the
                % link information. 
                num_new_endpoint = num_new_endpoint + 1;
%                 fprintf('Node %d only connects to single link. Convert it to be an endpoint\n', connected_node_label);
                [input_graph, new_link_info] = fun_graph_pruning_convert_node_to_endpoint(input_graph, connected_node_label);
                new_endpoint_ind(num_new_endpoint) = new_link_info.ep_ind;
                link_label = [];
            else
                link_label = [];
            end
        end
    end
end
new_endpoint_ind = new_endpoint_ind(1:num_new_endpoint);
%% Re-compute the labels
input_graph = fun_graph_relabel(input_graph);
end