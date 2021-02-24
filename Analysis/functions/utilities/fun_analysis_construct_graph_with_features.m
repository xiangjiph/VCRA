function region_graph = fun_analysis_construct_graph_with_features(region_graph, ...
    cube_link_feature_cell, cube_node_feature_cell)

fprintf('Collecting region links and nodes features\n');
region_feature.link = cube_link_feature_cell(cellfun(@(x) size(x, 1), cube_link_feature_cell) > 0);
region_feature.link = cat(1, region_feature.link{:});
region_graph.link.features = fun_analysis_get_link_features_from_local_stat_table(region_graph, region_feature.link);

region_feature.node = cube_node_feature_cell(cellfun(@(x) size(x, 1), cube_node_feature_cell) > 0);
region_feature.node = cat(1, region_feature.node{:});
region_graph.node.features = fun_analysis_get_node_features_from_local_stat_table(region_graph, region_feature.node);
end