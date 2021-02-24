vessel_graph = load('test_graph.mat');
max_capillary_r = 3.5; % um
internal_offset = 16; % um
%% Compute basic network features
vessel_graph.link.features = fun_analysis_get_link_features_by_label(vessel_graph);
vessel_graph.node.features = fun_analysis_get_node_cc_features_by_label(vessel_graph);
%% Initialization
current_graph = vessel_graph;
num_remove_link_per_step = 100;
num_cc_size_recorded = 10;
num_removal_steps = 50;
graph_cc_str_cell = cell(num_removal_steps, 1);
%%
for iter_step = 1 : num_removal_steps
    tmp_tic = tic;
    % Remove all the endpoints
    current_graph = fun_graph_pruning_internal_short_hairs(current_graph, inf, 0);
    tmp_link_graph_label = fun_analysis_get_link_graph_cc_label(current_graph);
    
    tmp_sub = cellfun(@(x) fun_ind2sub(current_graph.num.mask_size, x), current_graph.link.cc_ind, 'UniformOutput', false);
    tmp_link_length =  cellfun(@fun_graph_sub_to_length, tmp_sub);
    %% Record graph connected component data
    graph_cc_str = struct;
    graph_cc_str.link_length = tmp_link_length;
    graph_cc_str.total_link_length = sum(tmp_link_length);
    graph_cc_str.num_link = numel(tmp_link_graph_label);
    graph_cc_str.cc_link_label_cell = fun_bin_data_to_idx_list(tmp_link_graph_label);    
    graph_cc_str.num_cc = numel(graph_cc_str.cc_link_label_cell);
    graph_cc_str.num_link_per_cc = cellfun(@numel, graph_cc_str.cc_link_label_cell);
    [graph_cc_str.num_link_per_cc, tmp_sort_ind]= sort(graph_cc_str.num_link_per_cc, 'descend');
    graph_cc_str.cc_link_label_cell = graph_cc_str.cc_link_label_cell(tmp_sort_ind);
    % Largest connected component
    graph_cc_str.num_link_largest_cc = graph_cc_str.num_link_per_cc(1);
    graph_cc_str.largest_cc_ratio = graph_cc_str.num_link_largest_cc / graph_cc_str.num_link;
    [tmp_large_cc_num_links, tmp_large_cc_ratio] = deal(nan(1, num_cc_size_recorded));
    tmp_record_ind = 1 : min(num_cc_size_recorded, graph_cc_str.num_cc);
    tmp_large_cc_num_links(tmp_record_ind) = graph_cc_str.num_link_per_cc(tmp_record_ind);
    tmp_large_cc_ratio(tmp_record_ind) = tmp_large_cc_num_links(tmp_record_ind) ./ graph_cc_str.num_link;
    
    graph_cc_str.num_link_large_cc = tmp_large_cc_num_links;
    graph_cc_str.large_cc_ratio = tmp_large_cc_ratio;
    % Length
    graph_cc_str.cc_link_length_cell = cellfun(@(x) tmp_link_length(x), ...
        graph_cc_str.cc_link_label_cell, 'UniformOutput', false);
    graph_cc_str.cc_link_total_length = cellfun(@sum, graph_cc_str.cc_link_length_cell);
    graph_cc_str.largest_cc_total_length = graph_cc_str.cc_link_total_length(1);
    graph_cc_str.largest_cc_length_ratio = graph_cc_str.largest_cc_total_length / ...
        graph_cc_str.total_link_length;
    graph_cc_str_cell{iter_step} = graph_cc_str;
    %% Randomly remove links, as well as newly created endpoints
    tmp_rm_link_label = randsample(graph_cc_str.num_link, num_remove_link_per_step, false);
    % Do not worry about the endpoints here. The endpoint will be removed
    % by grpah pruning at the beginning of the iteration
    current_graph = fun_graph_delete_internal_links(current_graph, tmp_rm_link_label);
    fprintf('Finish link removal perturbation step %d. Elapsed time is %f seconds.\n', ...
        iter_step, toc(tmp_tic));
end
%% Analyze result
largest_cc_size = cellfun(@(x) x.num_link_largest_cc, ...
    graph_cc_str_cell);
largest_cc_length = cellfun(@(x) x.largest_cc_total_length, ...
    graph_cc_str_cell);
largest_cc_length_ratio = cellfun(@(x) x.largest_cc_length_ratio, ...
    graph_cc_str_cell);
%% 
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
yyaxis(ax_hdl, 'left');
plot(ax_hdl, largest_cc_length);
ax_hdl.YLabel.String = 'Total length in the largest cc(\mum)';
ax_hdl.XLabel.String = 'Iteration';
yyaxis(ax_hdl, 'right');
plot(ax_hdl, largest_cc_length_ratio);
ax_hdl.YLabel.String = 'Largest cc length ratio';

