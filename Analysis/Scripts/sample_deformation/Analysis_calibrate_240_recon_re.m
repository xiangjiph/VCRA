set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack_list = {'ML_2018_08_15', 'ML20190124', 'ML20200201'};
num_stack = numel(stack_list);
calibration_data = DataManager.load_data(DataManager.fp_metadata_file('Vessel_radius_calibration', ...
    'DK202004_merged', sprintf('Vessel_radius_calibration_DK202004_merged_mean_calibration_curve.mat')));
calibration_interpolation = calibration_data{end}.voxel.spline_itp;
%%
grid_version = '240_cube';
ori_recon_data_version = '240_cube_re';
output_recon_data_version = '240_cube_rec';
num_core = 12;
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
            tmp_cube_data = DataManager.load_block_skl(...
                tmp_grid_info.dataset_name, ...
                tmp_grid_info.stack, ori_recon_data_version, ...
                tmp_cube_grid_sub(1), tmp_cube_grid_sub(2), tmp_cube_grid_sub(3));
        catch ME
            fprintf('Unable to load data (label: %d)\n', iter_bbox);
            fprintf('Error message: %s\n', ME.message);
            continue;
        end
        try 
            tmp_cube_data = fun_analysis_sc_skeleton_radius(tmp_cube_data, calibration_interpolation);
            tmp_cube_data.grid_name = output_recon_data_version;
            
            DataManager.write_block_skl_file(tmp_cube_data, ...
                tmp_grid_info.dataset_name, ...
                tmp_grid_info.stack, output_recon_data_version,...
                tmp_cube_grid_sub(1), tmp_cube_grid_sub(2), tmp_cube_grid_sub(3));
        catch ME
            fprintf('Unable to process data (label: %d)\n', iter_bbox);
            rethrow(ME)
        end
        parfor_wait_message.Send;
    end    
    parfor_wait_message.Destroy;
    fprintf('Finish correcting scale of radius estimation in %s. Elapsed time is %f seconds.\n', ...
        stack, toc(tic_start));
end
pool_obj = gcp('nocreate');
delete(pool_obj);

