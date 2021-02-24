function vessel_graph = fun_analysis_get_whole_brain_graph(dataset_name, stack, grid_info, skl_version)

persistent DataManager
if isempty(DataManager)
    DataManager = FileManager;
end
%% Parallel load skeleton
wb_skel_str = fun_analysis_get_whole_brain_grid_skel(dataset_name, stack, grid_info, skl_version);
%%
whole_brain_skel_ind = cat(1, wb_skel_str.ind{:});
whole_brain_skel_r = cat(1, wb_skel_str.r{:});
[whole_brain_skel_ind_unique, unique_idx, ~] = unique(whole_brain_skel_ind);
whole_brain_skel_r_unique = whole_brain_skel_r(unique_idx);
% Takes 1400 seconds for conversion
tic
vessel_graph = fun_skeleton_to_graph(whole_brain_skel_ind_unique, grid_info.data_size);
toc
skel_r_sparse = sparse(whole_brain_skel_ind_unique, ones(size(whole_brain_skel_ind_unique)), double(whole_brain_skel_r_unique), vessel_graph.num.block_voxel, 1);
vessel_graph = fun_graph_add_radius(vessel_graph, skel_r_sparse, 0.25);
end