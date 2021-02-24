function [input_graph, undeleted_info]= fun_graph_delete_dim_short_link(input_graph, link_label_to_remove)
% fun_graph_delete_dim_short_link removes the links specified by the link
% label. If delecting the link result in creating a endpoint, this link
% will be skipped. 
%% Remove 
th_hair_length = 0;
num_link_to_remove = numel(link_label_to_remove);
new_endpoint_ind = zeros(ceil(num_link_to_remove/5),1);
num_new_endpoint = 0;
undeleted_info = fun_initialized_structure_array_with_fieldname_list({'change_lengthQ', ...
    'will_create_endpointQ', 'change_length_cc_ind', 'will_create_endpoint_cc_ind'});
if num_link_to_remove == 0
%     disp('No link needs to be removed');
    return;
else
%     fprintf('Need to remove %d links\n', num_link_to_remove);
end

input_graph.tmp_link_ind_old_to_new = 0 : input_graph.link.num_cc;
link_to_be_deleted_length = input_graph.link.num_voxel_per_cc(link_label_to_remove);
link_to_remove_cc = input_graph.link.cc_ind(link_label_to_remove);

undeleted_info.change_lengthQ = false(num_link_to_remove, 1);
undeleted_info.will_create_endpointQ = false(num_link_to_remove, 1);

for iter_link = 1 : num_link_to_remove
    link_label = link_label_to_remove(iter_link);  
    link_label = input_graph.tmp_link_ind_old_to_new(link_label + 1);
    if input_graph.link.num_voxel_per_cc(link_label) ~= link_to_be_deleted_length(iter_link)
%         fprintf('The length of the link %d to be removed has changed. Skip this link\n',...
%             link_label);
        undeleted_info.change_lengthQ(iter_link) = true;
        continue;
    end
    while ~isempty(link_label) && (link_label ~= 0)
%         connected_node_label = input_graph.link.connected_node_label(link_label,:);
        % Remove the link:
        [input_graph, updated_node_info] = fun_graph_refine_remove_single_link(input_graph, link_label, false, false);
        if isempty(updated_node_info.label)
            % This link is not deleted
            link_label = [];
            undeleted_info.will_create_endpointQ(iter_link) = true;
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
undeleted_info.change_length_cc_ind = link_to_remove_cc(undeleted_info.change_lengthQ);
undeleted_info.will_create_endpoint_cc_ind = link_to_remove_cc(undeleted_info.will_create_endpointQ);
%% Re-compute the labels
input_graph = fun_graph_relabel(input_graph);
end