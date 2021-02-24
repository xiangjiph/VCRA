clc;clear;
DataManager = FileManager;

script_root_folder = DataManager.SCRIPT_PATH;
dataset_name = 'WholeBrain';
stack = 'ML20200201';
grid_version = '240_cube';
grid_info = DataManager.load_grid(dataset_name, stack, grid_version);
%% Set task name
task_name = sprintf('Segmentation_%s', datestr(now, 'YYYYmmdd'));
task_folder = DataManager.fp_task_folder(dataset_name, stack, task_name);
if ~isfolder(task_folder)
    mkdir(task_folder);
end
task_str = struct;
task_str.DataManager = DataManager;
task_str.dataset_name = dataset_name;
task_str.stack = stack;
task_str.grid_version = grid_version;
task_str.grid_c_info = grid_info;
task_str.task_option = [];
task_str.task_function_name = 'fun_task_segmentation';
task_str.task_name = task_name;
%% For 1um voxel size
voxel_length_um = round(mean(grid_info.voxel_size_um));
seg_parameters = struct;
seg_parameters.voxel_length_um = round(mean(grid_info.voxel_size_um));
seg_parameters.data_type = grid_info.data_type;
seg_parameters.rod_filter_radius_um = 1;
seg_parameters.rod_filter_length_um = round(6*seg_parameters.rod_filter_radius_um + 1);
seg_parameters.rod_filter_num_omega = 6;
seg_parameters.vesselness.DoG_scale_list_um = [0.5, 1, 2]./seg_parameters.voxel_length_um;
seg_parameters.vesselness_th = 0.1;
seg_parameters.adp_th_scale_1_um = 8;
seg_parameters.adp_th_scale_2_um = 16;
seg_parameters.morp_min_cc_size = 27;
seg_parameters.max_pool_size = 8;
seg_parameters.min_bg_std = 250;
seg_parameters.grid_info = grid_info;

task_str.opt = seg_parameters;
%% Task distributed information
dist_info_cell = cell(2, 0);
dist_info_cell(:, end+1) = {'MACHINE 1', 1};
dist_info_cell(:, end+1) = {'MACHINE 1', 1};
dist_info_cell(:, end+1) = {'MACHINE 1', 2};
dist_info_cell(:, end+1) = {'MACHINE 1', 2};
dist_info_cell(:, end+1) = {'MACHINE 1', 2};
dist_info_cell(:, end+1) = {'MACHINE 1', 2};
dist_info_cell(:, end + 1) = {'MACHINE 2', 1};
dist_info_cell(:, end + 1) = {'MACHINE 2', 1};
dist_info_cell(:, end + 1) = {'MACHINE 2', 1};
dist_info_cell(:, end + 1) = {'MACHINE 2', 1};

task_str.overwrite_Q = false;
num_processor = size(dist_info_cell, 2);
%% Generate task structure automatically
num_subtask = num_processor;
bbox_list = 1 : grid_info.num_valid_cube;
num_total_bbox = numel(bbox_list);
avg_num_bbox_per_subtask = ceil(num_total_bbox/num_subtask);
begin_idx = 1 : avg_num_bbox_per_subtask : num_total_bbox;
end_idx = min(num_total_bbox, begin_idx + avg_num_bbox_per_subtask - 1);
%%
DataManager.sync_script_to_server;
for task_idx = 1 : num_subtask
    tmp_task_str = task_str;
    tmp_task_str.machine_name = dist_info_cell{1, task_idx};
    tmp_task_str.gpuDevice = dist_info_cell{2, task_idx};
    tmp_task_str.task_list = bbox_list(begin_idx(task_idx) : end_idx(task_idx));
    [tmp_cmd, tmp_task_str] = fun_task_get_task_str(tmp_task_str, task_idx, true);
    DataManager.run_command_on_machine(tmp_task_str.machine_name, tmp_task_str.task_cmd_bs);
    %     system(tmp_cmd)
end
%% Debug
% fun_task_segmentation(tmp_task_str);
