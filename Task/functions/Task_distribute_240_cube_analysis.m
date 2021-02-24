%% Input parameters for distributed work 
clc;clear;
set_env
DataManager = FileManager;
dataset_name = 'WholeBrain';
% stack = 'ML_2018_08_15';
stack = 'ML20190124';
% stack = 'ML20200201';
skl_grid_name = '240_cube_rec';
grid_c_version = '240_cube_combined_5_o_2';
grid_c = DataManager.load_grid(dataset_name, stack, grid_c_version);
output_mask_name = '240_cube_recon';
output_graph_name = '240_cube_combined_5_o_2';
if ~isfield(grid_c, 'internal_subgrid_label_array')
    grid_c = fun_grid_get_internal_subgrid(grid_c);
    DataManager.write_grid_info(grid_c, grid_c.dataset_name, grid_c.stack, grid_c.version);
end
task_function_name = 'fun_task_analysis_240_cube_reconstruction';
%% Basic information
task_str = struct;
task_str.dataset_name = dataset_name;
task_str.stack = stack;
task_str.DataManager = FileManager;
task_str.grid_c_name = grid_c_version;
task_str.grid_c_info = DataManager.load_grid(task_str.dataset_name, ...
    task_str.stack, task_str.grid_c_name);
task_str.task_function_name = task_function_name;
%% Task specific parameters
    %% For analyzing 240 cube
    opt_graph_feature = struct;
    computeQ = struct;
    computeQ.basic_link_feature = true;
    computeQ.basic_node_feature = true;
    computeQ.basic_reconstruction = true;
    computeQ.node_path_to_nearest_neighbor = true;
    computeQ.link_shortest_path = true;
    computeQ.dist_tissue_2_vessel = true;
    computeQ.dimension = false;
    computeQ.max_z_proj = true;
    computeQ.capillary_branching_order = true;
    computeQ.capillary_artery_order = NaN; % Placeholder for furture development 
    computeQ.capillary_vein_order = NaN; % Placeholder for furture development
    
    computeQ.noncapillary_DT = true; % Pick the vessels are are not capillaries, analyze the distance transform properties. 
    computeQ.sgl_seg_rm_ptb = true;
    computeQ.dist_to_brain_surface = true; % Compute the distance between the vessels and the surface of the brain. 
    
    computeQ.ep2ep_anisotropy = true;
    computeQ.volume_weighted_ep2ep_anisotropy = true;
    computeQ.capillary_ep2ep_anisotropy = true;
    computeQ.capillary_volume_weighted_ep2ep_anisotropy = true;    
    
    opt_graph_feature.computeQ = computeQ;
    opt_graph_feature.vis_dim_fit = false;
    opt_graph_feature.dim_fit_cutoff_length = 25;
    opt_graph_feature.recon_max_error_rate = 0.1;
    opt_graph_feature.merge_neighbor_nodes_Q = false;
    opt_graph_feature.max_merge_length = 5;
    opt_graph_feature.always_merge_max_length = 1;
    opt_graph_feature.DT_scale_factor = 0.5;
    opt_graph_feature.nonmax_win_size_um = 40;
    
    opt_graph_feature.capillary_max_radius = 3.5; % For saperating capillary vs large vessels for local statistics. 
    
    opt_graph_feature.brain_mask_fp = fullfile(DataManager.fp_mask_folder(dataset_name, stack, 'whole_brain_d16x_registration'), ...
        sprintf('%s_%s_d16x_registration_mask.nii.gz', dataset_name, stack));
    opt_graph_feature.brain_mask_downsample_rate = 16;
    assert(isfile(opt_graph_feature.brain_mask_fp), 'Brain mask file does not exist.');
    opt_graph_feature.brain_mask = niftiread(opt_graph_feature.brain_mask_fp) > 0;
    
    
    % Overwrite the existing features in the graph (if exist)
    opt_graph_feature.overwrite_computed_featureQ = false;
    
    task_str.skl_grid_name = skl_grid_name;
    task_str.subgrid_mask_name = output_mask_name;
    task_str.output_graph_name = output_graph_name;
    task_str.output_mask_name = output_mask_name;
    task_str.task_name = sprintf('Analysis_240_cube_reconstruction_%s', datestr(now, 'yyyymmdd'));
    task_str.graph_feature = opt_graph_feature;
%% Save default parameters for internal grid analysis
% DataManager.write_task_parameters(opt_graph_feature, 'default', 'internal_subgrid_analysis.xml');
%% Computation setting
dist_info_cell = cat(2, repelem({'MACHINE 1'}, 1, 4), repelem({'MACHINE 2'}, 1, 8));
task_str.overwrite_Q = false;
num_processor = size(dist_info_cell, 2);
random_shuffle_task_Q = true;
%% Distribute task
task_list = 1 : grid_c.num_valid_cube;
num_task = numel(task_list);
avg_num_task_per_process = ceil(num_task/num_processor);
begin_idx = 1 : avg_num_task_per_process : num_task;
end_idx = min(num_task, begin_idx + avg_num_task_per_process - 1);
assert(end_idx(end) == num_task);
if random_shuffle_task_Q
    if isrow(task_list)
        task_list = task_list.';
    elseif size(task_list, 1) < size(task_list, 2)
        warning('Number of row is less than number of column in the task list');
    end
    task_list = task_list(randperm(size(task_list,1)), :);
%     assert(numel(unique(task_list)) == grid_c.num_valid_cube);
end
%% Generate task structure and launch computation
DataManager.sync_script_to_server;
for iter_task = 1 : num_processor
    tmp_task_str = task_str;
    tmp_task_str.gpuDevice = [];
    tmp_task_str.machine_name = dist_info_cell{iter_task};
    tmp_task_str.task_list = task_list(begin_idx(iter_task) : end_idx(iter_task) , :);
    [tmp_cmd, tmp_task_str] = fun_task_get_task_str(tmp_task_str, iter_task, true);
    DataManager.run_command_on_machine(tmp_task_str.machine_name, tmp_task_str.task_cmd_bs);
    fprintf('Finish deploying task %d\n', iter_task);
end
fprintf('Finish distributing task.\n');
%%
% [tmp_cmd, tmp_task_str] = fun_task_get_task_str(task_str, 1, false);
% % tmp = fun_task_analysis_240_cube_reconstruction(tmp_task_str);
% % tmp = fun_task_analysis_240_cube_reconstruction('/data/Vessel/WholeBrain/ML_2018_08_15/task/Analysis_240_cube_reconstruction_20190702/Analysis_240_cube_reconstruction_20190702_processor_18.mat');
% task_str.overwrite_Q = true;
% tic
exit_code = fun_analysis_internal_subgrid(task_str.dataset_name, task_str.stack, ...
    task_str.skl_grid_name, 1947, task_str);
% toc