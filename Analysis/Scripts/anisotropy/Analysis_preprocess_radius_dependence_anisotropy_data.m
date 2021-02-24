set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack_list = {'ML_2018_08_15', 'ML20190124', 'ML20200201'};
num_stack = numel(stack_list);
stack_name_list = cellfun(@(x) strrep(x, '_', ''), stack_list, 'UniformOutput', false);
grid_version = '240_cube';
data_folder = 'Radius_dependence_anisotropy';
skel_version = '240_cube_rec';
num_core = 12;
%% Parameters for computation
weight_method = {'volume', 'length'};
%%
for iter_stack = 1 : num_stack
    stack = stack_list{iter_stack};
    grid_info = DataManager.load_grid(dataset_name, stack, grid_version);
    num_valid_cube = grid_info.num_valid_cube;
    %% Load data
    [wb_vw_data, wb_lw_data] = deal(cell(num_valid_cube, 1));
    load_data_tic = tic;
    pool_obj = gcp('nocreate');
    delete(pool_obj);
    parpool(num_core);
    parfor_wait_message = parfor_wait(num_valid_cube, ...
        'Waitbar', true);
    parfor iter_bbox = 1 : num_valid_cube
        try
            tmp_cube_data = DataManager.load_analysis_data_in_grid(data_folder, ...
                dataset_name, stack, skel_version, iter_bbox);
        catch ME
            fprintf('Unable to load cube %d. Error message: %s\n', iter_bbox, ME.message);
            continue;
        end
        try
            wb_lw_data{iter_bbox} = tmp_cube_data.length_weighted;
            wb_vw_data{iter_bbox} = tmp_cube_data.volume_weighted;
        catch ME
            fprintf('Unable to process cube %d. Error message:\n%s\n', iter_bbox, ME.message);
        end
        parfor_wait_message.Send;
    end
    pool_obj = gcp('nocreate');
    delete(pool_obj);
    parfor_wait_message.Destroy;
    fprintf('Finish loading data. Elapsed time is %f seconds\n', toc(load_data_tic));
    %% Process data
    % Find the nonempty structure
    nonempty_Q = ~cellfun(@isempty, wb_vw_data);
    radius_selection_min = cellfun(@(x) x.select_r_min, wb_vw_data(nonempty_Q), 'UniformOutput', false);
    radius_selection_min = unique(cat(1, radius_selection_min{:}), 'rows', 'stable');
    assert(size(radius_selection_min, 1) == 1, 'Multiple selection minimum vector');
    radius_selection_max = cellfun(@(x) x.select_r_max, wb_vw_data(nonempty_Q), 'UniformOutput', false);
    radius_selection_max = unique(cat(1, radius_selection_max{:}), 'row', 'stable');
    assert(size(radius_selection_max, 1) == 1, 'Multiple selection maximum vector');
    num_range = numel(radius_selection_max);
    assert(num_range == numel(radius_selection_min), 'Mismatch number of radius selection maximum and minimum');
    %%
    wb_data_str = struct;
    wb_data_str.datast_name = dataset_name;
    wb_data_str.stack = stack;
    wb_data_str.data_folder = data_folder;
    wb_data_str.skel_version = skel_version;
    wb_data_str.grid_version = grid_version;
    wb_data_str.radius_min = radius_selection_min;
    wb_data_str.radius_max = radius_selection_max;
    wb_data_str.num_radius_range = num_range;
    wb_data_str.filepath = DataManager.fp_analysis_data_file(dataset_name, stack, ...
        sprintf('%s_%s_%s_%s_stat_data.mat', dataset_name, stack, skel_version, data_folder));
    wb_data_str.volume_weighted = fun_analysis_collect_weighted_anisotropy_data_in_cell_array(...
        wb_vw_data, num_range);
    wb_data_str.length_weighted = fun_analysis_collect_weighted_anisotropy_data_in_cell_array(...
        wb_lw_data, num_range);
    DataManager.write_data(wb_data_str.filepath, wb_data_str);
    fprintf('Finish processing radius dependence anisotropy data for %s.\n', stack);
end