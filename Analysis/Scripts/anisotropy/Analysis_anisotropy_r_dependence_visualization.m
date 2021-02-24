set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
% stack = 'ML20190124';
% stack = 'ML20200201';
stack = 'ML_2018_08_15';
mask_name = '240_cube_recon_sc';
image_grid_name = '240_cube';
grid_c_version = '240_cube_combined_5_o_2';
grid_c_info = DataManager.load_grid(dataset_name, stack, grid_c_version);
grid_info = grid_c_info.grid_ori;
save_folder_name = 'whole_brain_stat_sc';
linear_scaling_factor = 1.0521;
registration_version = 'Allen_2017_25um_landmark.mat';
registration_str = DataManager.load_registration_data(dataset_name,...
    stack, registration_version);
skel_version = '240_cube_rec';
anisotropy_data_name = 'Radius_dependence_anisotropy';
vis_weight_method = 'volume_weighted';
%% Load data 
wb_max_proj_str = DataManager.load_analysis_data(dataset_name, stack, sprintf('%s_%s_%s_max_proj_data.mat', ...
    dataset_name, stack, mask_name), save_folder_name);
wb_240_cube_stat_str = DataManager.load_analysis_data(dataset_name, stack, sprintf('%s_%s_%s_240_cube_stat_data.mat', ...
    dataset_name, stack, mask_name), save_folder_name);
rda_data = DataManager.load_data(DataManager.fp_analysis_data_file(...
    dataset_name, stack, sprintf('%s_%s_%s_%s_stat_data.mat', dataset_name, ...
    stack, skel_version, anisotropy_data_name)));

wb_vis_anisotropy_data = rda_data.(vis_weight_method);
fprintf('Finish loading all the data\n');
%% Generate regional contour
% old regional ID
vis_region_id = {184, 453, 500, 1057, 677, 669, 247, 31, 972, 44, 714, 95, 254, 22, 541, 922, 895, ... % isocortex
    507, 151, 159, 589, 814, 961, 619, 631, 788, 566, ... % olfactory 
    375, 726, 982, 19, 822, ... % 
    703, ...
    477, 803, ...
    549, 1097, ...
    339, 323, 348, ...
    771, 354, ...
    528, 519, ... 
    776};  % Corpus callsoum
vis_region_id = cellfun(@(x) full(registration_str.map_oldID_to_newID(x)), vis_region_id, 'UniformOutput', false);
vis_region_name = cellfun(@(x) registration_str.label_2_name(x), ...
    vis_region_id, 'UniformOutput', false);
num_vis_region = numel(vis_region_id);

combined_regional_mask = zeros(registration_str.array_size, 'uint16');
for iter_region = 1 : num_vis_region
    % Get regional CC 
    tmp_region_id = vis_region_id{iter_region};
    [vis_region_cc]= fun_registration_get_region_cc(registration_str, tmp_region_id);
    if numel(tmp_region_id) > 1
        tmp_region_id = tmp_region_id(1);
    end
    combined_regional_mask(vis_region_cc.PixelIdxList{1}) = tmp_region_id;
    combined_regional_mask(vis_region_cc.PixelIdxList{2}) = tmp_region_id;
end

registration_vis_str = registration_str;
    registration_vis_str.combined_region_id = vis_region_id;
registration_vis_str.combined_region_name = vis_region_name;
registration_vis_str.combined_region_abbrv = cellfun(@(x) registration_str.label_2_acronym(x), ...
    vis_region_id, 'UniformOutput', false);
%% Generate cell array for visualization - anisotropy
% Exclude cube that has less than 25% of the volume outside the brain 
is_vis_cube_Q = (wb_240_cube_stat_str.cube_in_brain_mask_ratio >= 0.1);
vis_radius_range_idx = 10;
vis_radius_min = rda_data.radius_min(vis_radius_range_idx);
vis_radius_max = rda_data.radius_max(vis_radius_range_idx);
vis_filename_template = sprintf('%s_%%s_r_%.2f_to_%.2f_um', vis_weight_method, ...
    vis_radius_min, vis_radius_max);
anisotropy_plot_cell = cell(4, 0);

anisotropy_plot_cell(:, end+1) = {wb_vis_anisotropy_data.fa(:, vis_radius_range_idx), ...
    sprintf(vis_filename_template, 'FA'), 'Fractional Anisotropy', wb_vis_anisotropy_data.svd_max_vec(:, :, vis_radius_range_idx)};

anisotropy_plot_cell(:, end+1) = {wb_vis_anisotropy_data.fa_z(:, vis_radius_range_idx), ...
    sprintf(vis_filename_template, 'FAz'), 'FA_z', wb_vis_anisotropy_data.svd_max_vec(:, :, vis_radius_range_idx)};
%% Generate videos to visualize the basic statistics on three planes
for vis_dir = 1 : 3
    proj_im_str = fun_vis_get_recon_max_proj_str(wb_max_proj_str, vis_dir);
    switch vis_dir
        case 1
            registration_vis_str.combined_regional_mask = permute(combined_regional_mask, ...
                [2, 3, 1]);
        case 2
            registration_vis_str.combined_regional_mask = permute(combined_regional_mask, ...
                [1, 3, 2]);
        case 3
            registration_vis_str.combined_regional_mask = combined_regional_mask;
    end
        
    for iter_cell = 1 : size(anisotropy_plot_cell, 2)
        save_layer_imageQ = false;
        tmp_plot_cell = anisotropy_plot_cell(:, iter_cell);
        tmp_plot_cell{1}(~is_vis_cube_Q) = nan;
        fun_vis_generate_local_anisotropy_w_ori_video_w_boundary(tmp_plot_cell{:}, ...
            grid_info, proj_im_str, registration_vis_str, 'vec_proj', save_layer_imageQ, 5, 95);
    end
end
%% 
tmp_debug_array = nan(grid_info.grid_size);
tmp_debug_array(grid_info.bbox_grid_ind_list) = wb_vis_anisotropy_data.fa_z(:, vis_radius_range_idx);
implay(isfinite(tmp_debug_array))

debug_grid_sub = [28 35 40];
debug_cube_label = grid_info.bbox_grid_label_array(debug_grid_sub(1), ...
    debug_grid_sub(2), debug_grid_sub(3));
debug_anisotropy_data = DataManager.load_analysis_data_in_grid(anisotropy_data_name, ...
    dataset_name, stack, skel_version, debug_cube_label);