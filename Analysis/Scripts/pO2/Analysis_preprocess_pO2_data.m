set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack_list = {'ML_2018_08_15', 'ML20190124', 'ML20200201'};
% stack_list = {'ML20200201'};
stack_name_list = cellfun(@(x) strrep(x, '_', ''), stack_list, 'UniformOutput', false);
grid_version = '240_cube';
pO2_data_folder = 'pO2_SA_li_sc';

skel_version = '240_cube_rec';
local_extrema_window_size = [8 : 4 : 44, 50 : 10 : 70, 80 : 20 : 120];

%%
for iter_stack = 1 : numel(stack_list)
    stack = stack_list{iter_stack};
    tmp_grid_info = DataManager.load_grid(dataset_name, stack, grid_version);
    tic_start = tic;
    wb_pO2_stat = fun_analysis_load_whole_brain_pO2_dt_data(tmp_grid_info, pO2_data_folder);
    wb_pO2_data = fun_analysis_preprocess_pO2_data(wb_pO2_stat, numel(local_extrema_window_size));
    wb_pO2_data.dataset_name = dataset_name;
    wb_pO2_data.stack = stack;
    wb_pO2_data.skel_version = skel_version;
    wb_pO2_data.pO2_data_folder = pO2_data_folder;
    wb_pO2_data.local_extrema_window_size = local_extrema_window_size;
    wb_pO2_data.filepath = fullfile(DataManager.fp_analysis_data_folder(dataset_name, stack), ...
        sprintf('%s_%s_%s_%s_stat_data.mat', dataset_name, stack, skel_version, pO2_data_folder));
    DataManager.write_data(wb_pO2_data.filepath, wb_pO2_data);
    fprintf('Finish preprocessing whole brain pO2 simulation result in %s. Elapsed time is %f seconds.\n', ...
        stack, toc(tic_start));
end
