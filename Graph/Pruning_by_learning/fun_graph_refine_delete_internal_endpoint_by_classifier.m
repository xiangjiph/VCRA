function [input_graph, refine_info_str ] = fun_graph_refine_delete_internal_endpoint_by_classifier(input_graph, ...
    vessel_image, vessel_mask_dt, classifier)

% Parameters
internal_offset = 16;
pruning_max_length = 2;
%% Remove short hairs
input_graph = fun_graph_pruning_internal_short_hairs(input_graph, pruning_max_length, internal_offset);
%% Get internal points
int_link_0 = fun_graph_get_free_link(input_graph, internal_offset, true);
tmp_ep1_features = fun_graph_get_link_w_1_ep_features(int_link_0.ep1.link_cc_ind, vessel_image, vessel_mask_dt, true);
% Classification
tmp_ep1_link_removed_Q = classifier.classifier.predict(table2array(tmp_ep1_features(:, classifier.used_feature_idx)));
%
tmp_ep1_link_removed_label = int_link_0.ep1.link_label(tmp_ep1_link_removed_Q);
% Update graph 
input_graph = fun_graph_pruning_by_link_label(input_graph, tmp_ep1_link_removed_label);

refine_info_str.num_removed = numel(tmp_ep1_link_removed_label);
refine_info_str.kept_cc_ind = int_link_0.ep1.link_cc_ind(~tmp_ep1_link_removed_Q);
end