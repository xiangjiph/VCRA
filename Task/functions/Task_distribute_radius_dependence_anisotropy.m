%% Input parameters for distributed work 
clc;clear;
set_env
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
skl_grid_name = '240_cube_rec';
grid_version = '240_cube';

task_function_name = 'fun_task_radius_dependence_anisotropy';
save_folder_name = 'Radius_dependence_anisotropy';
% task_name = sprintf('%s_%s', save_folder_name, datestr(now, 'YYYYmmDD'));
task_name = sprintf('%s_20200603', save_folder_name);
%% Basic information
task_str = struct;
task_str.dataset_name = dataset_name;
task_str.stack = stack;
task_str.DataManager = DataManager;
task_str.grid_name = grid_version;
task_str.grid_info = DataManager.load_grid(task_str.dataset_name, ...
    task_str.stack, task_str.grid_name);
task_str.task_function_name = task_function_name;
task_str.task_name = task_name;
task_str.skl_grid_name = skl_grid_name;
%% Task specific parameters
    %% Simulation setting
    task_specific_opt = struct;
    radius_selection_max = [2 : 0.1 : 2.5, 2.75 : 0.25 : 3.5, 4, 5, 6, 8, 10, 200];
    num_radius_selection_range = numel(radius_selection_max);
    radius_selection_min = zeros(size(radius_selection_max));
    radius_selection_range = cat(1, radius_selection_min, radius_selection_max);
    radius_selection_range = cat(2, radius_selection_range, [3.5; 200]);
    task_specific_opt.radius_selection_range = radius_selection_range;
    task_specific_opt.weight_method = {'volume', 'length'};
    task_specific_opt.save_folder_name = save_folder_name;
    
    task_str.opt = task_specific_opt;
%% Task distributed infomation
dist_info_cell = cat(2, repelem({'MACHINE 1'}, 1, 2), repelem({'MACHINE 2'}, 1, 8));
task_str.overwrite_Q = false;
num_processor = size(dist_info_cell, 2);
random_shuffle_task_Q = true;
%% Determine if the combined grid contains valid internal cubes
% There's no point to solve the Poisson equation for regions on the surface
%% Distribute task
task_list = 1 : task_str.grid_info.num_valid_cube;
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
    tmp_task_str.machine_name = dist_info_cell{1, iter_task};
    tmp_task_str.task_list = task_list(begin_idx(iter_task) : end_idx(iter_task) , :);
    [tmp_cmd, tmp_task_str] = fun_task_get_task_str(tmp_task_str, iter_task, true);
    DataManager.run_command_on_machine(tmp_task_str.machine_name, tmp_task_str.task_cmd_bs);
    fprintf('Finish deploying task %d / %d.\n', iter_task, num_processor);
%     system(tmp_cmd)
end
%%
% fun_task_radius_dependence_anisotropy(tmp_task_str);