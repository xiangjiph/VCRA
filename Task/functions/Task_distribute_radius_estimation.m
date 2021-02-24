set_env
clc;clear;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML20200201';
grid_version = '240_cube';
load_skel_version = sprintf('%s_auto', grid_version);
write_skel_version = sprintf('%s_re', grid_version);
% Generate grid for radius estimation
grid_info = DataManager.load_grid(dataset_name, stack, grid_version);
%% Get the relative path
if ~isfield(grid_info.octree, 'relative_filepath')
    grid_info.octree.relative_filepath = cell(grid_info.octree.grid_size);
    current_root_path = DataManager.ROOTPATH;
    for iter_cell = 1 : grid_info.octree.num_block
        tmp_ind = grid_info.octree.grid_pos_ind(iter_cell);
        tmp_str = grid_info.octree.filepath{tmp_ind};
        if ~isempty(tmp_str)
            tmp_str = erase(tmp_str, sprintf('%s%s', current_root_path, '/'));
            grid_info.octree.relative_filepath{tmp_ind} = tmp_str;
        end
    end
    fprintf('Finish generating relative filepath for the rendering images\n');
    DataManager.write_grid_info(grid_info, dataset_name, stack, grid_version);
end
%%
re_grid = fun_grid_generate_combined_grid(grid_info, [5,5,1], 0, true);
task_function_name = 'fun_task_radius_estimation';
%% Basic information
task_str = struct;
task_str.dataset_name = dataset_name;
task_str.stack = stack;
task_str.DataManager = FileManager;
task_str.grid_version = grid_version;
task_str.task_function_name = task_function_name;
task_str.task_name = sprintf('Radius_estimation_%s', datestr(now, 'yyyymmdd'));
%% Task specific parameters
% Name
re_opt = struct;
re_opt.dataset_name = dataset_name;
re_opt.stack = stack;
re_opt.load_skel_version = load_skel_version;
re_opt.write_skel_version = write_skel_version;
% Data
re_opt.grid_info = DataManager.load_grid(dataset_name, stack, grid_version);
% Running
re_opt.parallel_load_Q = false;
% Radius estimation related
radius_est_grid_name = sprintf('%s_c%d_%d_%d_o_%d', grid_version, [5,5,1], 0);
re_grid = DataManager.load_grid(dataset_name, stack, radius_est_grid_name);
re_opt.grid_radius_estimation = re_grid;

psf_est_int = DataManager.load_data(DataManager.fp_metadata_file(dataset_name, ...
    stack, 'Estimated_PSF_edge_intensity'));
est_psf_int = psf_est_int.n_min_edge_int_interpolation;
re_opt.psf_edge_int_interpolation = est_psf_int;

re_opt.max_r_to_est_um = 15;
re_opt.num_ori_vec_vxl_half = 2;
re_opt.num_iteration = 6;
% maximum value of (estimated radius - original radius). Radius estimation
% greater than this value will be rejected.
re_opt.max_ext_to_ori_r_um = 2;
% Applied median filter to remove line-shift artefact
re_opt.render_data_medfilt3_Q = true;

task_str.radius_estimation_opt = re_opt;
%% Task distributed information
dist_info_cell = cat(2, repelem({'MACHINE 2'}, 1, 6), repelem({'MACHINE 1'}, 1, 12));
task_str.overwrite_Q = false;
num_processor = size(dist_info_cell, 2);
random_shuffle_task_Q = true;
%% Distribute task
task_list = 1 : re_grid.num_valid_cube;
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
    task_list = task_list(randperm(numel(task_list)));
    assert(numel(unique(task_list)) == re_grid.num_valid_cube);
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
%% Debug
% fun_task_radius_estimation(tmp_task_str);
% grid_c_list_ind = 2647;
% exit_code = fun_radius_estimation_in_combined_grid(re_opt, grid_c_list_ind);
