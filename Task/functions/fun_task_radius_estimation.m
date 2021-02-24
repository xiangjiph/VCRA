function exit_code = fun_task_radius_estimation(task_str)

if ~isa(task_str, 'struct')
    if isfile(task_str)
        task_str = load(task_str);
    else
        error('The input should be either task structure or filepath to the task structure');
    end
end
%% Phrase parameters
dataset_name = task_str.dataset_name;
stack = task_str.stack;

task_list = task_str.task_list;
num_task = numel(task_list);
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
%% Processing
diary(task_str.log_file_path);
fprintf('%s Start processing\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
record_t_start = tic;
for iter_grid_c = 1 : num_task
    diary('on');
    tmp_tic = tic;
    task_label = task_list(iter_grid_c);
    try 
        tmp_complete_log_fp = fullfile(task_str.task_record_folder, sprintf('%s_%s_%s_job_%d.txt', ...
            task_str.dataset_name, task_str.stack, task_str.task_name, task_label));
        if isfile(tmp_complete_log_fp) && ~task_str.overwrite_Q
            fprintf('Task has been completed previously. Do not overwrite.\n');
            continue;
        end
        exit_code = fun_radius_estimation_in_combined_grid(task_str.radius_estimation_opt, task_label);
        
        exit_info_str = sprintf('Finish processing 240-cube %d. Elapsed time is %f seconds. Exit code %d\n', task_label, toc(tmp_tic), exit_code);
        fprintf(exit_info_str);
        % Write a log file to the task folder
        system_write(exit_info_str, tmp_complete_log_fp, 'text');
    catch ME
        system_write(sprintf('Fail to process combined grid %d', task_label), ...
            task_str.error_file_path, 'text');
        system_write(sprintf('Error message: %s', getReport(ME, 'extended', 'hyperlinks', 'off')), ...
            task_str.error_file_path, 'text');        
    end
    diary('off');
end
diary('on');
fprintf('%s Finish task. Elapsed time is %f seconds\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), toc(record_t_start));
system(sprintf('rm %s', task_str.filepath));
fprintf('Finish deleting task structure file.\n');
diary('off');
exit_code = 0;
end