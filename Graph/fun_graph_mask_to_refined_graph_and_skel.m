function exit_code = fun_graph_mask_to_refined_graph_and_skel(dataset_name, stack, grid_c_version, combined_grid_xyz_label, opt)


%% Parameters
save_combined_grid_Q = opt.save_combined_grid_Q;
overwrite_existing_graphQ = opt.overwrite_existing_graphQ;
overwrite_skl_Q = opt.overwrite_skl_Q;
%% Initialization
DataManager = FileManager;
if isfield(opt, 'grid_c')
    grid_c = opt.grid_c;
else
    grid_c = DataManager.load_grid(dataset_name, stack, grid_c_version);
end
%%
data_info = fun_graph_get_grid_c_block_info(grid_c, combined_grid_xyz_label);
if isfield(opt, 'output_graph_name') && ~isempty(opt.output_graph_name)
    data_info.grid_name = opt.output_graph_name;
end
if isfield(opt, 'output_skel_name') && ~isempty(opt.output_skel_name)
    data_info.sub_grid_version = opt.output_skel_name;
end
% Check if the graph exist
graph_fp = DataManager.fp_graph_in_block_file(dataset_name, stack, data_info.grid_name, data_info.idx_1, ...
    data_info.idx_2, data_info.layer);
if isfile(graph_fp) && ~overwrite_existing_graphQ
    fprintf('Have already processed %s. Skip\n', graph_fp);
    exit_code = 1;
    return;
end
%% Load data
combined_grid_mmxx_grid = data_info.combined_grid_mmxx_grid;
vessel_mask = DataManager.load_blocks_files('mask', dataset_name, stack, grid_c.grid_ori.version, ...
    combined_grid_mmxx_grid(1):combined_grid_mmxx_grid(4), combined_grid_mmxx_grid(2):combined_grid_mmxx_grid(5), ...
    combined_grid_mmxx_grid(3):combined_grid_mmxx_grid(6), 'logical');

vessel_image = DataManager.load_blocks_files('image', dataset_name, stack, grid_c.grid_ori.version, ...
    combined_grid_mmxx_grid(1):combined_grid_mmxx_grid(4), combined_grid_mmxx_grid(2):combined_grid_mmxx_grid(5), ...
    combined_grid_mmxx_grid(3):combined_grid_mmxx_grid(6), 'uint16');
%% Convert skeleton to graph
[vessel_graph, vessel_mask_dt] = fun_graph_mask_to_graph(vessel_mask, vessel_image, opt.mask2graph);
%% Correction and refinement
% try
ori_vessel_skel_ind = cat(1, vessel_graph.link.pos_ind, vessel_graph.node.pos_ind);
for iter_cycle = 1 : opt.auto_refine_cycle
    [vessel_graph, opt.record] = fun_graph_refine_by_classifier(vessel_graph, vessel_image, vessel_mask_dt, opt);
end
refined_vessel_skel_ind = cat(1, vessel_graph.link.pos_ind, vessel_graph.node.pos_ind);
opt.record.deleted_skel_ind = setdiff(ori_vessel_skel_ind, refined_vessel_skel_ind);
vessel_graph.record = opt.record;
% catch ME
%     fprintf('Fail to refine combined block %d. Error message: %s\n', combined_grid_xyz_label, ME.identifier);
% end
%% Save graph
vessel_graph.info = data_info;
vessel_graph = fun_graph_add_radius(vessel_graph, vessel_mask_dt, true);
if save_combined_grid_Q
    DataManager.write_graph_in_block(vessel_graph, dataset_name, stack, data_info.grid_name, ...
        data_info.idx_1, data_info.idx_2, data_info.layer);
end
exit_code = fun_graph_write_subblock_skeleton(vessel_graph, overwrite_skl_Q);
end