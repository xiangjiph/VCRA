function output_table = fun_analysis_get_link_features_from_local_stat_table(region_graph, link_feature)

rm_feature_list = {'mid_global_ind', 'ep1_global_ind', 'ep2_global_ind', ...
    'link_com', 'ep1_sub', 'ep2_sub', 'ep1_direction_vec', 'ep2_direction_vec', ...
    'cc_sub_pca1_vec', 'cc_sub_pca2_vec', 'cc_sub_pca3_vec', ...
    'cc_sub_cov_eig_val', 'mid_sub', 'mid_global_sub', 'dt_max', 'dt_min', ...
    'shortest_loop_node_label', 'shortest_loop_link_label', 'ep2ep_angle_azimuth_deg', ...
    'ep2ep_angle_elevation_deg'};
%%
ep_pair_ind = cat(2, link_feature.ep1_global_ind, link_feature.ep2_global_ind);
ep_pair_ind = sort(ep_pair_ind, 2, 'ascend');
[unique_ep_pair_ind, unique_ep_list_ind, ~] = unique(ep_pair_ind, 'row', 'stable');
is_valid_ep_pair_Q = all(~isnan(unique_ep_pair_ind), 2);
unique_ep_pair_ind = unique_ep_pair_ind(is_valid_ep_pair_Q, :);
unique_ep_list_ind = unique_ep_list_ind(is_valid_ep_pair_Q);
% Link label of the endpoint pairs in the region_graph
ep_pair_link_label = full(region_graph.link.map_ind_2_label(unique_ep_pair_ind));
% 
is_in_the_region_graph_Q = all(ep_pair_link_label > 0, 2);
% Select the consistent ones
is_in_the_same_link_Q = (ep_pair_link_label(:, 1) == ep_pair_link_label(:, 2)) & ...
    is_in_the_region_graph_Q;
% Check if the endpoint pairs that are not in the region_graph are all on
% the edge of the ROI. The two endpoints should be out of the ROI ? 

% Check if the endpoint pairs that are in different link of region_graph
% are all on the edge of the ROI? These links go out of the ROI and come
% back again. 

% Determine which row of features to add to the region_graph.link.features
ep_pair_link_label_same = ep_pair_link_label(is_in_the_same_link_Q, 1); % region_graph link label
link_label_same_ep_pair_ind = unique_ep_pair_ind(is_in_the_same_link_Q, :);
link_label_same_ep_list_ind = unique_ep_list_ind(is_in_the_same_link_Q); % Feature list ind
[unique_ep_pair_same_link_list_ind, ep_pair_link_label_same_u] = fun_bin_data_to_idx_list(ep_pair_link_label_same);
% Determine which feature to use for the given endpoint pairs that are in
% the same link in the region_graph 
unique_ep_pair_same_link_list_ind_vec = nan(size(ep_pair_link_label_same_u));
for iter_link = 1 : numel(unique_ep_pair_same_link_list_ind)
    tmp_list_ind = unique_ep_pair_same_link_list_ind{iter_link};
    tmp_num_list_elem = numel(tmp_list_ind);
    if tmp_num_list_elem == 1
        unique_ep_pair_same_link_list_ind_vec(iter_link) = tmp_list_ind;
    elseif tmp_num_list_elem > 1
        % Check if the endpoint pairs are the endpoint of the link
        % connected component in the region_graph 
        tmp_ep_pair_list = link_label_same_ep_pair_ind(tmp_list_ind, :);
        tmp_link_label = ep_pair_link_label_same_u(iter_link);
        tmp_link_cc = region_graph.link.cc_ind{tmp_link_label};
        tmp_match_Q = false;
        for iter_pair = 1 : tmp_num_list_elem
            tmp_ep_pair = tmp_ep_pair_list(iter_pair, :);
            tmp_ep1_ind = find(tmp_link_cc(1) == tmp_ep_pair);
            if ~isempty(tmp_ep1_ind)
                if ~isscalar(tmp_ep1_ind)
                    % Is it a correct way to do? 
                    tmp_ep1_ind = tmp_ep1_ind(1);
                end
                switch tmp_ep1_ind
                    case 1
                        if tmp_link_cc(end) == tmp_ep_pair(2)
                            % Found
                            unique_ep_pair_same_link_list_ind_vec(iter_link) = tmp_list_ind(iter_pair);
                            tmp_match_Q = true;
                            break;
                        end
                    case 2
                        if tmp_link_cc(end) == tmp_ep_pair(1)
                            % Found
                            unique_ep_pair_same_link_list_ind_vec(iter_link) = tmp_list_ind(iter_pair);
                            tmp_match_Q = true;
                            break;
                        end
                    otherwise
                        error('The endpoint pairs are not the endpoints of the link cc in the region_graph');
                end
            end
        end
        if ~tmp_match_Q
            assert(all(any(bsxfun(@eq, tmp_link_cc', tmp_ep_pair_list(:)), 2)), 'Not all the endpoint pairs are in the link cc');
        end
    elseif tmp_num_list_elem < 1
        error('Empty cell');
    end
end
% If is not nan, the corresponding link features can be copied to the
% region_graph link feature. 
valid_list_ind_Q_1 = ~isnan(unique_ep_pair_same_link_list_ind_vec);
valid_feature_list_ind_1 = link_label_same_ep_list_ind(unique_ep_pair_same_link_list_ind_vec(valid_list_ind_Q_1));
valid_link_label_1 = ep_pair_link_label_same_u(valid_list_ind_Q_1);

% Add the link features to the link 
% Why some link in the graph does not have features? 
feature_name = link_feature.Properties.VariableNames;
feature_name = setdiff(feature_name, rm_feature_list);
num_feature = numel(feature_name);
output_table = struct;
for iter_feature = 1 : num_feature
    tmp_feature_name = feature_name{iter_feature};
    tmp_feature_value = link_feature(valid_feature_list_ind_1, tmp_feature_name);
    tmp_feature_value = table2array(tmp_feature_value);
    if iscell(tmp_feature_value)
        tmp_value = cell(region_graph.link.num_cc, 1);
        tmp_value(valid_link_label_1, :) = tmp_feature_value;
    elseif isvector(tmp_feature_value)
        tmp_value = nan(region_graph.link.num_cc, 1);
        tmp_value(valid_link_label_1) = tmp_feature_value;
    elseif ismatrix(tmp_feature_value)
        tmp_value = nan(region_graph.link.num_cc, size(tmp_feature_value, 2));
        tmp_value(valid_link_label_1, :) = tmp_feature_value;        
    end
    output_table.(tmp_feature_name) = tmp_value;
end
output_table = struct2table(output_table);
%% Check if the calculation is correct
% tmp_graph_link_length = nan(region_graph.link.num_cc, 1);
% for iter_link = 1 : region_graph.link.num_cc
%     tmp_graph_link_length(iter_link) = fun_graph_sub_to_length(fun_ind2sub(region_graph.num.mask_size, ...
%         region_graph.link.cc_ind{iter_link}));
% end
% tmp_is_valid_Q = ~isnan(region_graph.link.features.length);
% tmp_consistent_Q = (output_table.length(tmp_is_valid_Q) - tmp_graph_link_length(tmp_is_valid_Q)) < 1e-2;
% if (all(tmp_consistent_Q))
%     fprintf('Check pass. All the link in the region graph is not longer than the link length computed in the local analysis.\n');
% %     fprintf('Fraction of links with identical length: %f\n', ...
% %         nnz(output_table.length(tmp_is_valid_Q) == tmp_graph_link_length(tmp_is_valid_Q))/numel(tmp_consistent_Q));
% else
%     fig_hdl = figure;
%     ax_hdl = axes(fig_hdl);
%     scatter(ax_hdl, tmp_graph_link_length, output_table.length)
%     grid(ax_hdl, 'on');
%     ax_hdl.XLabel.String = 'Link length computed in the region graph/\mum';
%     ax_hdl.YLabel.String = 'Link length loaded from regional analysis/\mum';
%     ax_hdl.DataAspectRatio = [1,1,1];
%     error('The link length computed in the region_graph is longer than the length computed in the local analysis\n');
% end
% All the datapoint in the scatter plot should be on or under y = x curve
% up to numerical error
end

%% Clear up 
