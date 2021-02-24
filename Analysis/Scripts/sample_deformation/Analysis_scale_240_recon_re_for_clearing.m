% For consistency and convenience, the voxel size of the rendered image
% were scale up directly by the linear shrinkage ratio due to tissue
% clearing. However, evidences suggest that vessel radius is largely
% unaffected by the tissue clearing. Therefore, the estiamted radius should
% be scaled back before going into graph analysis part. 
set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack_list = {'ML20200201'};
num_stack = numel(stack_list);
clearing_scaling_factor = 1.0926;
radius_scaling_factor = 1 / clearing_scaling_factor;
%%
grid_version = '240_cube';
ori_recon_data_version = '240_cube_re';
output_recon_data_version = '240_cube_re';
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
            tmp_cube_data = fun_analysis_sc_skeleton_radius(tmp_cube_data, radius_scaling_factor);
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

