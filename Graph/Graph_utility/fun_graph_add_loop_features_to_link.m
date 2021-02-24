function link_feature = fun_graph_add_loop_features_to_link(link_feature, tmp_loop_str)


% Initialize the feature
num_link = size(link_feature, 1);
link_feature.shortest_loop_lenght = nan(num_link, 1);
link_feature.shortest_loop_geodesic_lenght = nan(num_link, 1);
link_feature.shortest_loop_link_label = cell(num_link, 1);
link_feature.shortest_loop_node_label = cell(num_link, 1);

% Add features according to the link label
link_feature.shortest_loop_lenght(tmp_loop_str.link_label) = tmp_loop_str.loop_length;
link_feature.shortest_loop_geodesic_lenght(tmp_loop_str.link_label) = tmp_loop_str.loop_geodesic_length;
link_feature.shortest_loop_link_label(tmp_loop_str.link_label) = tmp_loop_str.loop_link_label;
link_feature.shortest_loop_node_label(tmp_loop_str.link_label) = tmp_loop_str.loop_node_label;
%% Generative derivative features
link_feature.avg_edge_length_in_shortest_loop = link_feature.shortest_loop_lenght ...
    ./ link_feature.shortest_loop_geodesic_lenght;

link_feature.length_ratio_in_shortest_loop = link_feature.length ./ ...
    link_feature.shortest_loop_lenght;
end