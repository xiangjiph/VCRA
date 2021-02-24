set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack_list = {'ML_2018_08_15', 'ML20190124', 'ML20200201'};
% stack_list = {'ML20200201'};
linear_scaling_factor = 1.0521;
radius_scaling_factor = 1;
% clearing_scaling_factor = 1.0926;
% radius_scaling_factor = linear_scaling_factor / clearing_scaling_factor;
%% Correct pO2 data and save to new folder
pO2_grid_version = '240_cube';
ori_pO2_data_folder = 'pO2_SA_li';
output_pO2_data_folder = sprintf('%s_sc', ori_pO2_data_folder);
num_core = 10;
for iter_stack = 1 : numel(stack_list)
    stack = stack_list{iter_stack};
    tmp_grid_info = DataManager.load_grid(dataset_name, stack, pO2_grid_version);
    tic_start = tic;    
    pool_obj = gcp('nocreate');
    delete(pool_obj);
    parpool(num_core);
    parfor_wait_message = parfor_wait(tmp_grid_info.num_valid_cube, ...
        'Waitbar', true);
    parfor iter_bbox = 1 : tmp_grid_info.num_valid_cube
        % Load data
        try 
            tmp_cube_data = DataManager.load_analysis_data_in_grid(...
                ori_pO2_data_folder, tmp_grid_info.dataset_name, ...
                tmp_grid_info.stack, tmp_grid_info.version, iter_bbox);
        catch ME
            fprintf('Unable to load data (label: %d)\n', iter_bbox);
            fprintf('Error message: %s\n', ME.message);
            continue;
        end
        try 
            tmp_cube_data = fun_analysis_sc_pO2_result(tmp_cube_data, linear_scaling_factor, radius_scaling_factor);
            DataManager.write_analysis_data_in_grid(tmp_cube_data, ...
                output_pO2_data_folder, tmp_grid_info.dataset_name, ...
                tmp_grid_info.stack, tmp_grid_info.version, iter_bbox);
        catch ME
            fprintf('Unable to process data (label: %d)\n', iter_bbox);
            rethrow(ME)
        end
        parfor_wait_message.Send;
    end    
    parfor_wait_message.Destroy;
    fprintf('Finish correcting scale of whole brain pO2 simulation result in %s. Elapsed time is %f seconds.\n', ...
        stack, toc(tic_start));
end
pool_obj = gcp('nocreate');
delete(pool_obj);
%% Correct reconstruction data
% The radius is not scaled 
grid_version = '240_cube';
ori_recon_version = '240_cube_recon';
output_recon_version = sprintf('%s_sc', ori_recon_version);
num_core = 10;
for iter_stack = 1 : numel(stack_list)
    stack = stack_list{iter_stack};
    tmp_grid_info = DataManager.load_grid(dataset_name, stack, grid_version);
    tic_start = tic;    
    pool_obj = gcp('nocreate');
    delete(pool_obj);
    parpool(num_core);
    parfor_wait_message = parfor_wait(tmp_grid_info.num_valid_cube, ...
        'Waitbar', true);
    parfor iter_bbox = 1 : tmp_grid_info.num_valid_cube
        tmp_cube_grid_sub = tmp_grid_info.bbox_grid_sub_list(iter_bbox, :);
        % Load data
        try 
            tmp_cube_data = DataManager.load_block_mask(...
                tmp_grid_info.dataset_name, tmp_grid_info.stack, ori_recon_version, ...
                tmp_cube_grid_sub(1), tmp_cube_grid_sub(2), tmp_cube_grid_sub(3));
        catch ME
            fprintf('Unable to load data (label: %d)\n', iter_bbox);
            fprintf('Error message: %s\n', ME.message);
            continue;
        end
        try 
            tmp_cube_data = fun_analysis_sc_reconstruction_result(tmp_cube_data, linear_scaling_factor, radius_scaling_factor);
            tmp_cube_data.version = output_recon_version;
            DataManager.write_block_mask(tmp_cube_data, ...
                tmp_grid_info.dataset_name, tmp_grid_info.stack, ...
                output_recon_version, tmp_cube_grid_sub(1), tmp_cube_grid_sub(2), tmp_cube_grid_sub(3));
        catch ME
            fprintf('Unable to process data (label: %d)\n', iter_bbox);
            rethrow(ME)
        end
        parfor_wait_message.Send;
    end    
    parfor_wait_message.Destroy;
    fprintf('Finish correcting scale of whole brain pO2 simulation result in %s. Elapsed time is %f seconds.\n', ...
        stack, toc(tic_start));
end
pool_obj = gcp('nocreate');
delete(pool_obj);
%% Correct the length density
grid_version = '240_cube';
skl_grid_name = '240_cube_rec';
correct_field_name = {'all', 'cap', 'all_w_cep', 'all_w_eph', 'all_w_ep', ...
    'cap_w_cep', 'cap_w_eph', 'cap_w_ep'};
for iter_stack = 2 : numel(stack_list)
    stack = stack_list{iter_stack};
    cube_length_density_fp = fullfile(DataManager.fp_analysis_data_folder(dataset_name, ...
        stack), sprintf('%s_%s_%s_cube_length_density.mat', ...
        dataset_name, stack, skl_grid_name));
    cube_len_den_str = DataManager.load_data(cube_length_density_fp);
    for iter_fn = 1 : numel(correct_field_name)
        tmp_fn = correct_field_name{iter_fn};
        cube_len_den_str.(tmp_fn) = cube_len_den_str.(tmp_fn) ./ (linear_scaling_factor^2);
    end
    cube_len_den_str.filepath = fullfile(DataManager.fp_analysis_data_folder(dataset_name, ...
        stack), sprintf('%s_%s_%s_cube_length_density_sc.mat', ...
        dataset_name, stack, skl_grid_name));
    DataManager.write_data(cube_len_den_str.filepath, cube_len_den_str);
end
