function fun_task_segmentation_v1(task_str)
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
if ~isfield(task_str, 'gpuDevice')
    gpuDevice(2);
else
    gpuDevice(task_str.gpuDevice);
end
if ~isfield(task_str, 'fun_handle')
    task_fun = @fun_mouselight_segmentation_1um_cube;
else
    task_fun = task_str.fun_handle;
end

diary(task_str.log_file_path);
record_t_start = tic;
DataManager = task_str.DataManager;
fprintf('%s\tSegmentating task start\n', datestr(now));
for tmp_block_idx = 1 : task_str.task_num_bbox
    diary('on');
    grid_idx_1 = task_str.task_bbox_list(tmp_block_idx,1);
    grid_idx_2 = task_str.task_bbox_list(tmp_block_idx,2);
    grid_layer = task_str.task_bbox_list(tmp_block_idx,3);
    fprintf('Processing cube (%d, %d, %d)\n', grid_idx_1, grid_idx_2, grid_layer);
    try
        block_data = DataManager.load_block_data(task_str.dataset_name, task_str.stack, task_str.grid_version, grid_idx_1, grid_idx_2, grid_layer);
        mask_str = struct;
        [vessel_mask, mask_str.record] = task_fun(block_data, task_str.seg_parameters);
        mask_str.ind = uint32(find(vessel_mask));
        mask_str.block_size = size(vessel_mask);
        DataManager.write_block_mask(mask_str, task_str.dataset_name, task_str.stack, task_str.grid_version, ...
            grid_idx_1, grid_idx_2, grid_layer);
    catch ME
        system_write(sprintf('Fail to process cube (%d, %d, %d)\n', grid_idx_1, grid_idx_2, grid_layer), ...
            task_str.error_file_path, 'text');
        system_write(sprintf('Error message: %s', ME.identifier), ...
            task_str.error_file_path, 'text');               
    end
    diary('off');
end
diary('on');
fprintf('Finish task. Elapsed time is %f seconds\n', toc(record_t_start));
diary('off');
end