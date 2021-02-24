%% Input parameters for distributed work 
% To do list: 
% 1. Merge the neighboring node. Debug. 
clc;clear;
set_env
DataManager = FileManager;
dataset_name = 'WholeBrain';
% stack = 'ML20200201';
stack = 'ML20190124';
% stack = 'ML_2018_08_15';
skl_grid_name = '240_cube_rec';
grid_c_version = '240_cube_combined_5_o_2';
grid_c = DataManager.load_grid(dataset_name, stack, grid_c_version);
if ~isfield(grid_c, 'internal_subgrid_label_array')
    grid_c = fun_grid_get_internal_subgrid(grid_c);
    DataManager.write_grid_info(grid_c, grid_c.dataset_name, grid_c.stack, grid_c.version);
end
task_function_name = 'fun_task_simulation_pO2';
% task_name = sprintf('%s_%s', 'Simulation_pO2_SA', datestr(now, 'YYYYmmDD'));
task_name = sprintf('%s_%s', 'Simulation_pO2_SA', '20200726');
%% Basic information
task_str = struct;
task_str.dataset_name = dataset_name;
task_str.stack = stack;
task_str.DataManager = DataManager;
task_str.grid_c_name = grid_c_version;
task_str.grid_c_info = DataManager.load_grid(task_str.dataset_name, ...
    task_str.stack, task_str.grid_c_name);
task_str.task_function_name = task_function_name;
task_str.overwrite_Q = false;
%% Task specific parameters
    %% Simulation setting
    task_specific_opt = struct;
    task_specific_opt.recon_max_error_rate = 0.1; % Reconstruction 
    task_specific_opt.recon_pad_half_length = 80; % Before downsampling
    task_specific_opt.local_extrema_window_size = [8 : 4 : 44, 50 : 10 : 70, 80 : 20 : 120]; % Before downsampling 
    task_specific_opt.recon_downsample_rate = 2; % Downsample the reconstruction to reduce computation load
    task_specific_opt.min_recon_vessel_radius_um = 0; % The minimum radius of the vessel used for reconstruction. Not sure if 2 is a more appropriate number. 
    % Poisson equation 
    task_specific_opt.inhomogeneous_term = 1; % Also before downsampling
    % Krogh model
    task_specific_opt.krogh_vessel_r_um = 2;
    task_specific_opt.krogh_est_corr_coeff = 0.5;
    task_specific_opt.mask_itp_method = 'linear';
    task_specific_opt.use_gpu_Q = true;
    task_specific_opt.use_scratch_Q = true;
    task_specific_opt.save_folder_name = 'pO2_SA_li';
    
    task_str.skl_grid_name = skl_grid_name;
    task_str.task_name = task_name;
    task_str.opt = task_specific_opt;
%% Task distributed infomation
dist_info_cell = cell(2, 0);
dist_info_cell(:, end + 1) = {'MACHINE 1', 1};
dist_info_cell(:, end + 1) = {'MACHINE 1', 1};
dist_info_cell(:, end + 1) = {'MACHINE 1', 2};
dist_info_cell(:, end + 1) = {'MACHINE 1', 2};
% dist_info_cell(:, end + 1) = {'MACHINE 1', 2};
% dist_info_cell(:, end + 1) = {'MACHINE 1', 2};
% dist_info_cell(:, end + 1) = {'MACHINE 2', 1};
% dist_info_cell(:, end + 1) = {'MACHINE 2', 1};
% dist_info_cell(:, end + 1) = {'MACHINE 2', 1};
% dist_info_cell(:, end + 1) = {'MACHINE 2', 1};
% dist_info_cell(:, end + 1) = {'MACHINE 2', 1};
% dist_info_cell(:, end + 1) = {'MACHINE 2', 1};
% dist_info_cell(:, end + 1) = {'MACHINE 2', 1};
dist_info_cell(:, end + 1) = {'MACHINE 2', 1};
dist_info_cell(:, end + 1) = {'MACHINE 2', 1};
dist_info_cell(:, end + 1) = {'MACHINE 2', 2};
dist_info_cell(:, end + 1) = {'MACHINE 2', 2};
%% Computation setting
num_processor = size(dist_info_cell, 2);
random_shuffle_task_Q = true;
%% Determine if the combined grid contains valid internal cubes
% There's no point to solve the Poisson equation for regions on the surface
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
% Synchronize the data
if task_specific_opt.use_scratch_Q
    unique_computer_name = unique(dist_info_cell(1, :));
    for iter_machine = 1 : numel(unique_computer_name)
        tmp_machine_name = unique_computer_name{iter_machine};
        tmp_file_manager = FileManager(tmp_machine_name);
        tmp_source_folder_path = tmp_file_manager.fp_block_skl_folder(dataset_name, ...
            stack, task_str.skl_grid_name);
        tmp_target_folder_path = strrep(tmp_source_folder_path, tmp_file_manager.ROOT_PATH, ...
            tmp_file_manager.SCRATCH_ROOT_PATH);
        tmp_sync_str = sprintf('mkdir -p %s; rsync -rav %s/ %s', ...
            tmp_target_folder_path, tmp_source_folder_path, tmp_target_folder_path);
        DataManager.run_command_on_machine(tmp_machine_name, tmp_sync_str, false);
    end
end
%%
for iter_task = 1 : num_processor
    tmp_task_str = task_str;
    tmp_task_str.machine_name = dist_info_cell{1, iter_task};
    tmp_task_str.gpuDevice = dist_info_cell{2, iter_task};
    tmp_task_str.task_list = task_list(begin_idx(iter_task) : end_idx(iter_task) , :);
    [tmp_cmd, tmp_task_str] = fun_task_get_task_str(tmp_task_str, iter_task, true);
    DataManager.run_command_on_machine(tmp_task_str.machine_name, tmp_task_str.task_cmd_bs);
%     system(tmp_cmd)
    fprintf('Finish deploying job %d\n', iter_task)
end
%%
% exit_code = fun_task_simulation_pO2(tmp_task_str);
% exit_code = fun_simulation_OT_internal_subgrid_SA(task_str.grid_c_info, ...
%     task_str.skl_grid_name, 1459, task_str.opt);
%% Debug negative krogh fitting coefficient
% debug_240_cube_label = 9;
% debug_240_cube_grid_ind = grid_c.grid_ori.bbox_grid_ind_list(debug_240_cube_label);
% debug_240_cube_grid_c_label = grid_c.internal_subgrid_label_array(debug_240_cube_grid_ind);
% exit_code = fun_simulation_OT_internal_subgrid(task_str.grid_c_info, ...
%     task_str.skl_grid_name, debug_240_cube_grid_c_label, task_str.opt);
