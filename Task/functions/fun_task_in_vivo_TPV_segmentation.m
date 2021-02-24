function exit_code = fun_task_in_vivo_TPV_segmentation(task_str)
% function for distributing the segmentation task
% task_str = struct;
% task_str.dataset_name = dataset_name;
% task_str.stack = stack;
% task_str.grid_version = grid_version;
% task_str.seg_parameters = seg_parameters;
% task_str.output_log_file = output_log_file;
% task_str.DataManager = FileManager;
% 
% task_str.task_bbox_list = grid_info.bbox_grid_sub{1};
% task_str.task_num_bbox = numel(task_str.task_bbox_list(:,1));
if ~isa(task_str, 'struct')
    if isfile(task_str)
        task_str = load(task_str);
    else
        error('The input should be either task structure or filepath to the task structure');
    end
end
%% Parse parameters
dataset_name = task_str.dataset_name;
stack = task_str.stack;
grid_version = task_str.grid_version;
opt = task_str.opt;
task_list = task_str.task_list;
assert(isvector(task_list), 'task_list should be a vector');
num_task = numel(task_list);
% Set the gpuDevice
if isfield(task_str, 'gpuDevice')
    gpuDevice(task_str.gpuDevice);
else
    gpuDevice(2);
end

if isfield(task_str, 'DataManager')
    DataManager = task_str.DataManager;
else
    DataManager = FileManager;
end

grid_info = DataManager.load_grid(dataset_name, stack, grid_version);
%% Create folders
log_folder = fileparts(task_str.log_file_path);
if ~isfolder(log_folder)
    mkdir(log_folder);
end
error_folder = fileparts(task_str.error_file_path);
if ~isfolder(error_folder)
    mkdir(error_folder);
end
if ~isfolder(task_str.task_record_folder)
    mkdir(task_str.task_record_folder);
end
%%
diary(task_str.log_file_path);
record_t_start = tic;
fprintf('%s Start processing\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

for tmp_block_idx = 1 : num_task
    diary('on');
    task_label = task_list(tmp_block_idx);
    exit_code = -1;
    tmp_tic = tic;
    try
        tmp_complete_log_fp = fullfile(task_str.task_record_folder, sprintf('%s_%s_%s_job_%d.txt', ...
            task_str.dataset_name, task_str.stack, task_str.task_name, task_label));
        if isfile(tmp_complete_log_fp) && ~task_str.overwrite_Q
            fprintf('Task has been completed previously. Do not overwrite.\n');
            continue;
        end        
        % Unpack the function for performance. Packing the following script
        % into a function increase the computational time by 100%...
%         exit_code = task_fun(dataset_name, stack, grid_version, task_label, opt);
        grid_idx_1 = grid_info.bbox_grid_sub_list(task_label, 1);
        grid_idx_2 = grid_info.bbox_grid_sub_list(task_label, 2);
        grid_layer = grid_info.bbox_grid_sub_list(task_label, 3);

        mask_str = struct;
        mask_str.dataset_name = dataset_name;
        mask_str.stack = stack;
        mask_str.version = grid_version;
        mask_str.global_bbox_mmll = grid_info.bbox_xyz_mmll_list(task_label, :);
        mask_str.global_bbxx_mmxx = grid_info.bbox_xyz_mmxx_list(task_label, :);
        mask_str.global_block_size = grid_info.data_size;
        mask_str.grid_idx_1 = grid_idx_1;
        mask_str.grid_idx_2 = grid_idx_2;
        mask_str.grid_layer = grid_layer;
        % Load data
        block_data = DataManager.load_block_data(dataset_name, stack, grid_version, grid_idx_1, grid_idx_2, grid_layer);
        % Generate segmentation
        [vessel_mask, mask_str.record] = fun_in_vivo_TPV_segmentation_1um_cube(block_data, opt);
        % Convert mask to ind list
        assert(numel(vessel_mask) < intmax('uint32'));
        mask_str.ind = uint32(find(vessel_mask));
        mask_str.block_size = size(vessel_mask);
        % Record the max projection image for visualization
        mask_str.max_proj_1 = squeeze(max(block_data, [], 1));
        mask_str.max_proj_2 = squeeze(max(block_data, [], 2));
        mask_str.max_proj_3 = max(block_data, [], 3);
        % Save result
        DataManager.write_block_mask(mask_str, dataset_name, stack, grid_version, ...
            grid_idx_1, grid_idx_2, grid_layer);
        exit_code = 0;        
        
        exit_info_str = sprintf('Finish processing 240-cube %d. Elapsed time is %f seconds. Exit code %d\n', task_label, toc(tmp_tic), exit_code);
        fprintf(exit_info_str);
        % Write a log file to the task folder
        system_write(exit_info_str, tmp_complete_log_fp, 'text');
    catch ME
        system_write(sprintf('Fail to process task %d\n', task_label), ...
            task_str.error_file_path, 'text');
        system_write(sprintf('Error message: %s', getReport(ME, 'extended', 'hyperlinks', 'off')), ...
            task_str.error_file_path, 'text');               
    end
    diary('off');
end
diary('on');
fprintf('Finish task. Elapsed time is %f seconds\n', toc(record_t_start));
system(sprintf('rm %s', task_str.filepath));
fprintf('Finish deleting task structure file.\n');
diary('off');
exit_code = 0;
end