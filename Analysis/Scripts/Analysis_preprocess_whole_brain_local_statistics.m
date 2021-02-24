set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack_list = {'ML_2018_08_15', 'ML20190124', 'ML20200201'};
reconstruction_name = '240_cube_recon_sc';
grid_version = '240_cube';
save_folder_name = 'whole_brain_stat_sc';

% stack = 'ML_2018_08_15';
for iter_stack = 1 : numel(stack_list)
    tic_start = tic;
    stack = stack_list{iter_stack};
    exit_code = fun_analysis_collect_whole_brain_local_statistics(dataset_name, stack, grid_version, reconstruction_name, save_folder_name);
    fprintf('Finish preprocessing whole brain local statistics in %s. Elapsed time is %f seconds.\n', ...
        stack, toc(tic_start));
end
