clc;clear;close all;
vessel_graph = load('test_graph.mat');
vessel_graph_wo_ep = fun_graph_pruning_internal_short_hairs(vessel_graph, inf, 0);
%% Initialization
num_link = vessel_graph_wo_ep.link.num_cc;
link_rm_fraction_list = 0.05 : 0.05 : 0.90;
num_rm_fraction = numel(link_rm_fraction_list);
num_simulation = 20;
num_cc_size_recorded = 10;
graph_cc_str_cell = cell(num_simulation, num_rm_fraction);
%%
for iter_rm_fraction = 1 : num_rm_fraction
    tmp_rm_fraction = link_rm_fraction_list(iter_rm_fraction);
    tmp_num_rm_link_num = round(tmp_rm_fraction * vessel_graph_wo_ep.link.num_cc);
    for iter_simu = 1 : num_simulation
        tmp_tic = tic;
        %% Randomly remove links, as well as newly created endpoints
        tmp_rm_link_label = randsample(vessel_graph_wo_ep.link.num_cc, tmp_num_rm_link_num, false);
        tmp_kept_Q = true(num_link, 1);
        tmp_kept_Q(tmp_rm_link_label) = false;
        tmp_current_skel_ind = cat(1, vessel_graph_wo_ep.node.pos_ind, ...
            vessel_graph_wo_ep.link.cc_ind{tmp_kept_Q}, ...
            vessel_graph_wo_ep.isopoint.pos_ind);
        tmp_current_graph = fun_skeleton_to_graph(tmp_current_skel_ind, vessel_graph_wo_ep.num.mask_size);        
%         tmp_current_graph = fun_graph_delete_internal_links(vessel_graph_wo_ep, tmp_rm_link_label);
        tmp_current_graph = fun_graph_pruning_internal_short_hairs(tmp_current_graph, inf, 0);
%         Remove all the endpoints
        tmp_link_graph_label = fun_analysis_get_link_graph_cc_label(tmp_current_graph);
        %% Compute resulting graph properties
        % Compute the length of the remianing links 
        tmp_link_length = nan(tmp_current_graph.link.num_cc, 1);
        for iter_link = 1 : tmp_current_graph.link.num_cc
            tmp_link_length(iter_link) = fun_graph_sub_to_length(fun_ind2sub(...
                tmp_current_graph.num.mask_size, tmp_current_graph.link.cc_ind{iter_link}));
        end
        %% Compute distance transform properties
        vessel_recon_mask_label = fun_graph_to_reconstruction_mask(tmp_current_graph, true, 0.1);
        recon_prop = fun_analysis_reconstruction_space_properties(vessel_recon_mask_label, 0.5, false);
        %% Record graph connected component data
        graph_cc_str = struct;
        graph_cc_str.link_length = tmp_link_length;
        graph_cc_str.total_link_length = sum(tmp_link_length);
        graph_cc_str.num_link = numel(tmp_link_graph_label);
        graph_cc_str.cc_link_label_cell = fun_bin_data_to_idx_list(tmp_link_graph_label);
        graph_cc_str.num_cc = numel(graph_cc_str.cc_link_label_cell);
        graph_cc_str.num_link_per_cc = cellfun(@numel, graph_cc_str.cc_link_label_cell);
        [graph_cc_str.num_link_per_cc, tmp_sort_ind]= sort(graph_cc_str.num_link_per_cc, 'descend');
        graph_cc_str.total_num_link = sum(graph_cc_str.num_link_per_cc);
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
%%
        graph_cc_str_cell{iter_simu, iter_rm_fraction} = graph_cc_str;
        fprintf('Link removal fraction %f. Finish link removal perturbation step %d / %d. Elapsed time is %f seconds.\n', ...
            tmp_rm_fraction, iter_simu, num_simulation, toc(tmp_tic));
    end
end
%% Analyze result
largest_cc_size = cellfun(@(x) x.num_link_largest_cc, ...
    graph_cc_str_cell);
largest_cc_ratio = cellfun(@(x) x.largest_cc_ratio, ...
    graph_cc_str_cell);
largest_cc_length = cellfun(@(x) x.largest_cc_total_length, ...
    graph_cc_str_cell);
largest_cc_length_ratio = cellfun(@(x) x.largest_cc_length_ratio, ...
    graph_cc_str_cell);
%%
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2,1];
ax_hdl = subplot(1,2,1);
yyaxis(ax_hdl, 'left');
errorbar(ax_hdl, link_rm_fraction_list, mean(largest_cc_size, 1), std(largest_cc_size, 1));
ax_hdl.YLabel.String = 'Largest cc size';
ax_hdl.XLabel.String = 'Fraction of removed links';
yyaxis(ax_hdl, 'right');
errorbar(ax_hdl, link_rm_fraction_list, mean(largest_cc_ratio, 1), std(largest_cc_ratio, 1));
ax_hdl.YLabel.String = 'Largest cc ratio';
ax_hdl_2 = subplot(1,2,2);
yyaxis(ax_hdl_2, 'left');
errorbar(ax_hdl_2, link_rm_fraction_list, mean(largest_cc_length, 1), std(largest_cc_length, 1));
ax_hdl_2.YLabel.String = 'Total length in the largest cc(\mum)';
ax_hdl_2.XLabel.String = 'Fraction of removed links';
yyaxis(ax_hdl_2, 'right');
errorbar(ax_hdl_2, link_rm_fraction_list, mean(largest_cc_length_ratio, 1), std(largest_cc_length_ratio, 1));
ax_hdl_2.YLabel.String = 'Largest cc length ratio';
%%
dataset_name = vessel_graph.info.dataset_name;
stack = vessel_graph.info.stack;
grid_version = vessel_graph.info.grid_version;
% fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'percolation', ...
%     sprintf('%s_%s_%s_%d_giant_cc_vs_link_rm_fraction.png', dataset_name, stack, grid_version, vessel_graph.info.combined_grid_xyz_label));
fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'percolation', ...
    sprintf('%s_%s_%s_%d_giant_cc_vs_link_rm_fraction_wo_pruning.png', dataset_name, stack, grid_version, vessel_graph.info.combined_grid_xyz_label));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
