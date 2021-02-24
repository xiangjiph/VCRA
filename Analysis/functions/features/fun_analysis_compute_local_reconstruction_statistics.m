function cube_stat_str = fun_analysis_compute_local_reconstruction_statistics(wholebrain_stat_str)
% This implementation is too slow. Maybe because of broadcasting the entire
% cell feature array?
% Number of cores to run the computation in parallel:
num_core = 18;
%% Initialization
num_bbox = numel(wholebrain_stat_str.grid_ind);

[all_link_stat_cell, capillary_link_stat_cell, all_node_stat_cell] = deal(cell(num_bbox, 1));
%%
% Create cell array for parfor. Distributing wholebrain_stat_str is too
% expansive.
persistent uni_ori_isotropy_info
if isempty(uni_ori_isotropy_info)
    uni_ori_isotropy_info = load('./Metadata/uni_ori_isotropy.mat');
end

link_cell = wholebrain_stat_str.link_features;
node_cell = wholebrain_stat_str.node_features;
pool_obj = gcp('nocreate');
delete(pool_obj);
parpool(num_core);
parfor iter_bbox = 1 : num_bbox
    % for iter_bbox = 1 : num_bbox
    %     tmp_tic = tic;
    tmp_link_table = link_cell{iter_bbox};
    tmp_node_table = node_cell{iter_bbox};
    %%
    if ~isempty(tmp_link_table)
        % Statistics for all the vessels
        all_link_stat_cell{iter_bbox} = fun_analysis_get_graph_feature_basic_stat(tmp_link_table);
        % Statistics for capillary
        tmp_capillary_Q = ~tmp_link_table.is_large_vessel_Q;
        capillary_link_stat_cell{iter_bbox} = fun_analysis_get_graph_feature_basic_stat(tmp_link_table, tmp_capillary_Q);
    end
    % Statistics for all the nodes
    if ~isempty(tmp_node_table)
        all_node_stat_cell{iter_bbox} = fun_analysis_get_graph_feature_basic_stat(tmp_node_table);
    end
    %     fprintf('Finish processing cube %d. Elapsed time is %f seconds.\n', iter_bbox, toc(tmp_tic));
end
%%
pool_obj = gcp('nocreate');
delete(pool_obj);
cube_stat_str = struct;
cube_stat_str.all_link_stat = all_link_stat_cell;
cube_stat_str.capillary_stat = capillary_link_stat_cell;
cube_stat_str.all_node_stat = all_node_stat_cell;
end