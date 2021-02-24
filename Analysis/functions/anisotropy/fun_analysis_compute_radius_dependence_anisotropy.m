function exit_code = fun_analysis_compute_radius_dependence_anisotropy(grid_info, ...
    skel_version, cube_list_ind, opt)

persistent DataManager;
if isempty(DataManager)
    DataManager = FileManager;
end
%% Parse input
save_folder_name = opt.save_folder_name;
radius_selection_range = opt.radius_selection_range;
weight_method = opt.weight_method;

dataset_name = grid_info.dataset_name;
stack = grid_info.stack;
%%
% Load skeleton
tmp_grid_sub = grid_info.bbox_grid_sub_list(cube_list_ind, :);
try
    block_skl = DataManager.load_block_skl(dataset_name, stack, skel_version, ...
        tmp_grid_sub(1), tmp_grid_sub(2), tmp_grid_sub(3));
catch ME
    fprintf('Fail to load skeleton in grid %d %d %d. Terminate.\n', ...
        tmp_grid_sub);
    fprintf('Error message: %s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    exit_code = -1;
    return;
end
% Construct vessel graph
vessel_graph = fun_skeleton_to_graph(block_skl.ind, block_skl.block_size);
vessel_graph.radius = sparse(double(block_skl.ind), 1, double(block_skl.r), prod(block_skl.block_size), 1);
%% Initialize output structure
anisotropy_data = fun_getfields(block_skl, {'dataset_name', 'stack', 'grid_name', ...
    'idx_1', 'idx_2', 'layer'});
num_weight_method = numel(weight_method);
for iter_weight_method = 1 : num_weight_method
    tmp_weight_method = weight_method{iter_weight_method};
    anisotropy_all_vw_range = fun_analysis_get_anisotropy_stat_vs_radius_range(vessel_graph, ...
        radius_selection_range, tmp_weight_method);
    anisotropy_data.(sprintf('%s_weighted', tmp_weight_method)) = anisotropy_all_vw_range;
    % Post processing
end
%% Save data
DataManager.write_analysis_data_in_grid(anisotropy_data, ...
    save_folder_name, dataset_name, stack, skel_version, cube_list_ind);
exit_code = 0;
end