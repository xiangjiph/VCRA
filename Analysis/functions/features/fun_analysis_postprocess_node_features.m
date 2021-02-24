function node_features = fun_analysis_postprocess_node_features(node_features, delete_cell_Q)
% The features added here are all scale-invariant.
if nargin < 2
    delete_cell_Q = false;
end
% For the features stored in cell array
num_node = size(node_features, 1);
process_feature_name = {'conn_neig_dist_rank', 'link_ori_ep_cos', 'link_ep2ep_cos'};
save_feature_name_prefix = {'conn_node_drk', 'link_ori_agl', 'link_ep2ep_agl'};
for iter_feature = 1 : numel(process_feature_name)
    tmp_feature_name = process_feature_name{iter_feature};
    if ismember(tmp_feature_name, {'link_ori_ep_cos', 'link_ep2ep_cos'})
        convert_to_degree_Q = true;
    else
        convert_to_degree_Q = false;
    end
    % Initialize subfield
    tmp_name_min = sprintf('%s_min', save_feature_name_prefix{iter_feature});
    tmp_name_max = sprintf('%s_max', save_feature_name_prefix{iter_feature});
    tmp_name_med = sprintf('%s_med', save_feature_name_prefix{iter_feature});
    [node_features.(tmp_name_min), node_features.(tmp_name_med), ...
        node_features.(tmp_name_max)] = deal(nan(num_node, 1));
    
    for iter_node = 1 : num_node
        tmp_cell_data = node_features.(tmp_feature_name ){iter_node};
        tmp_cell_data = tmp_cell_data(isfinite(tmp_cell_data));
        if isempty(tmp_cell_data)
            continue;
        end
        if convert_to_degree_Q
            % Fix the numerical error in projection...
            tmp_cell_data(tmp_cell_data > 1) = 1;
            tmp_cell_data(tmp_cell_data < -1) = -1;
            tmp_cell_data = acosd(tmp_cell_data);
            assert(all(isreal(tmp_cell_data)));
        end
        tmp_cell_data = sort(tmp_cell_data , 'ascend');
        node_features.(tmp_name_max)(iter_node) = tmp_cell_data(end);
        node_features.(tmp_name_min)(iter_node) = tmp_cell_data(1);
        tmp_numel = numel(tmp_cell_data);
        if tmp_numel >= 2
            if mod(tmp_numel, 2) == 1
                node_features.(tmp_name_med)(iter_node) = tmp_cell_data((tmp_numel + 1) / 2);
            else
                node_features.(tmp_name_med)(iter_node) = 0.5 * (tmp_cell_data(tmp_numel/ 2) + ...
                    tmp_cell_data( tmp_numel/ 2 + 1));
            end
        end
    end
    if delete_cell_Q 
        node_features.(tmp_feature_name) = [];
    end
end
%% Delete unused field 
field_to_delete = {'nearest_node_label', 'center_of_mass', 'neighbor_node_label', 'distance_to_boundary'};
for iter_field = 1 : numel(field_to_delete)
    tmp_field_name = field_to_delete{iter_field};
    if ismember(tmp_field_name, node_features.Properties.VariableNames)
        node_features.(tmp_field_name) = [];
    end
end
end