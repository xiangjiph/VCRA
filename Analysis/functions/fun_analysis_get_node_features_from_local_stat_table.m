function [output_table, varargout]= fun_analysis_get_node_features_from_local_stat_table(region_graph, node_feature)


rm_field_name = {'nearest_node_label', 'center_of_mass', 'neighbor_node_label', ...
    'distance_to_boundary'};
%%
assert(isnumeric(node_feature.global_ind), 'The node global ind should be numeric');
[node_global_ind, unique_node_list_ind, ~] = unique(node_feature.global_ind, 'stable');

table_node_label = full(region_graph.node.map_ind_2_label(node_global_ind));
num_node_in_graph = region_graph.node.num_cc;

is_in_graphQ = table_node_label > 0;
fprintf('The number of nodes not in the region graph is %d\n', nnz(~is_in_graphQ));
% These nodes are probably on the boundary of the structure. 

% node_global_ind = node_global_ind(is_in_graphQ);
valid_node_list_ind = unique_node_list_ind(is_in_graphQ);
table_node_label = table_node_label(is_in_graphQ);

assert(numel(unique(table_node_label)) == numel(table_node_label), 'More than one nodes in the feautre table are mapped to the same node in the graph');

feature_name = node_feature.Properties.VariableNames;
feature_name = setdiff(feature_name, rm_field_name);
num_feature = numel(feature_name);
output_table = struct;
for iter_feature = 1 : num_feature
    tmp_feature_name = feature_name{iter_feature};
    tmp_feature_value = node_feature(valid_node_list_ind, tmp_feature_name);
    tmp_feature_value = table2array(tmp_feature_value);
    if iscell(tmp_feature_value)
        tmp_value = cell(num_node_in_graph, 1);
        tmp_value(table_node_label, :) = tmp_feature_value;
    elseif isvector(tmp_feature_value)
        tmp_value = nan(num_node_in_graph, 1);
        tmp_value(table_node_label) = tmp_feature_value;
    elseif ismatrix(tmp_feature_value)
        tmp_value = nan(num_node_in_graph, size(tmp_feature_value, 2));
        tmp_value(table_node_label, :) = tmp_feature_value;        
    end
    output_table.(tmp_feature_name) = tmp_value;
end

output_table = fun_struct2table(output_table);
%% Debug not in the graph node
if nargout > 1
    missing_node_global_ind = node_global_ind(~is_in_graphQ);
    missing_voxel_not_in_graph_Q = full(region_graph.radius(missing_node_global_ind)) == 0;
    missing_node_global_ind = missing_node_global_ind(missing_voxel_not_in_graph_Q);
    varargout{1} = missing_node_global_ind;
end
end