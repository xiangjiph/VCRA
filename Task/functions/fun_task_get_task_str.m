function [tmp_cmd, varargout] = fun_task_get_task_str(task_str, sub_task_ID, saveQ)
% fun_task_get_task_str generate the command line string for running the
% task by launching a new matlab session.

if nargin < 3
    saveQ = true;
end
persistent DataManager_local;
if isempty(DataManager_local)
    DataManager_local = FileManager;
end

if isfield(task_str, 'machine_name') && ~isempty(task_str.machine_name)
    DataManager_remote = FileManager(task_str.machine_name);
    script_root_folder = DataManager_remote.SERVER_SCRIPT_PATH;
else
    script_root_folder = DataManager_local.SERVER_SCRIPT_PATH;
end
cmd_format = 'cd %s; matlab -nodisplay -nosplash -r "addpath(genpath(''./'')); try, %s(''%s''); catch ME, %s quit; end; quit;" &';
error_handle_format = 'system_write(sprintf(''Error message: %%s\\n'', ME.message), ''%s'', ''text'');';
% error_handle_format = 'system_write(sprintf(''Error message: %%s\\n'', getReport(ME, ''extended'', ''hyperlinks'', ''off''), ''%s'', ''text'');';
%%
filename_ow_ext = sprintf('%s_processor_%d', task_str.task_name, sub_task_ID);
task_str.filename = sprintf('%s.mat', filename_ow_ext);
task_str.task_folder = DataManager_local.fp_task_folder(task_str.dataset_name, task_str.stack, task_str.task_name);
if ~isfolder(task_str.task_folder)
    mkdir(task_str.task_folder);
end
save_file_name = fullfile(task_str.task_folder, task_str.filename);
if isfield(task_str, 'machine_name')
     task_str.task_folder = strrep(task_str.task_folder, DataManager_local.ROOT_PATH, DataManager_remote.ROOT_PATH);
end
% Task setting filepath 
task_str.filepath = fullfile(task_str.task_folder, task_str.filename);
% Bash script filepath - directory in the remote computer
task_str.bs_filepath = fullfile(task_str.task_folder, sprintf('%s.sh', filename_ow_ext));
% Handle the task progress log
task_str.task_record_folder = fullfile(task_str.task_folder, 'completed');
% Handle error of the task
tmp_error_handle_str = sprintf(error_handle_format, fullfile(task_str.task_folder, 'task_error.txt'));
% Handle the standard output of the task 
tmp_log_folder = fullfile(task_str.task_folder, 'log');

task_str.log_file_name = sprintf('%s_processor_%d_log.txt', task_str.task_name, sub_task_ID);
task_str.log_file_path = fullfile(tmp_log_folder, task_str.log_file_name);
% Handle the error of the task function
tmp_error_folder = fullfile(task_str.task_folder, 'error_message');

task_str.error_file_name = sprintf('%s_processor_%d_error.txt', task_str.task_name, sub_task_ID);
task_str.error_file_path = fullfile(tmp_error_folder, task_str.error_file_name);
% Generate the command line string
tmp_cmd = sprintf(cmd_format, script_root_folder, task_str.task_function_name, task_str.filepath, tmp_error_handle_str);
task_str.task_cmd = tmp_cmd;
%%
if saveQ
    save(save_file_name, '-struct', 'task_str');
    % Write bash script
    % directory accessible to the local machine 
    [folder_name, file_name, ~] = fileparts(save_file_name); 
    bash_fp = fullfile(folder_name, sprintf('%s.sh', file_name));
    fid = fopen(bash_fp, 'w');
    fprintf(fid, '#!/bin/bash\n');
%     fprintf(fid, 'screen;\n');
    fprintf(fid, 'cd %s;\n', script_root_folder);
    fprintf(fid, 'nohup matlab -nodisplay -nosplash -r "addpath(genpath(''./'')); try, %s(''%s''); catch ME, %s quit; end; quit;" &\n', ...
        task_str.task_function_name, task_str.filepath, tmp_error_handle_str);
%     fprintf(fid, 'matlab -nodisplay -nosplash -r "addpath(genpath(''./'')); try, %s(''%s''); catch ME, %s quit; end; quit;" &\n', ...
%         task_str.task_function_name, task_str.filepath, tmp_error_handle_str);
%     fprintf(fid, 'screen -d;\n');
    fclose(fid);
    system(sprintf('chmod g+rx %s', bash_fp));
    task_str.task_cmd_bs = sprintf('source %s &', task_str.bs_filepath);
end
if nargout > 1
    varargout{1} = task_str;
end
end