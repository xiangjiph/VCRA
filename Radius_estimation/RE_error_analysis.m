set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML20190124';
mask_name = '240_cube_recon';
image_grid_name = '240_cube';
grid_info = DataManager.load_grid(dataset_name, stack, image_grid_name);
save_folder_name = 'whole_brain_stat';
%% Debug radius
stack_list = {'ML_2018_08_15', 'ML20190124'};
[wb_max_proj_str_cell, wb_240_cube_stat_str_cell, grid_info_cell] = deal(cell(size(stack_list)));
for iter_stack = 1 : numel(stack_list)
    tmp_stack = stack_list{iter_stack};
    grid_info_cell{iter_stack} = DataManager.load_grid(dataset_name, tmp_stack, image_grid_name);
    wb_max_proj_str_cell{iter_stack} = DataManager.load_analysis_data(dataset_name, tmp_stack, sprintf('%s_%s_%s_max_proj_data.mat', ...
        dataset_name, tmp_stack, mask_name), save_folder_name);
    wb_240_cube_stat_str_cell{iter_stack} = DataManager.load_analysis_data(dataset_name, tmp_stack, sprintf('%s_%s_%s_240_cube_stat_data.mat', ...
        dataset_name, tmp_stack, mask_name), save_folder_name);
end
%%
for iter_stack = 1 : numel(stack_list)
    proj_im_str = fun_vis_get_recon_max_proj_str(wb_max_proj_str_cell{iter_stack}, 3);
    fun_vis_generate_local_stat_video_fix_range(wb_240_cube_stat_str_cell{iter_stack}.link_cap_stat.median.dt_median,...
        'Capillary_median_radius', 'Radius(\mum)' , grid_info_cell{iter_stack}, proj_im_str, 0.75, 1.75, false);
end