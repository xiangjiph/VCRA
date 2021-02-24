set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
% stack = 'ML_2018_08_15';
stack = 'ML20200201';
mask_name = '240_cube_recon_sc';
image_grid_name = '240_cube';
grid_c_version = '240_cube_combined_5_o_2';
grid_c_info = DataManager.load_grid(dataset_name, stack, grid_c_version);
grid_info = grid_c_info.grid_ori;
save_folder_name = 'whole_brain_stat_sc';
linear_scaling_factor = 1.0521;
%% Load data 
wb_max_proj_str = DataManager.load_analysis_data(dataset_name, stack, sprintf('%s_%s_%s_max_proj_data.mat', ...
    dataset_name, stack, mask_name), save_folder_name);
wb_240_cube_stat_str = DataManager.load_analysis_data(dataset_name, stack, sprintf('%s_%s_%s_240_cube_stat_data.mat', ...
    dataset_name, stack, mask_name), save_folder_name);
% wb_pO2_stat_str = DataManager.load_analysis_data(dataset_name, stack, sprintf('%s_%s_%s_re_pO2_min_r_2um_sc_stat_data.mat', ...
%     dataset_name, stack, image_grid_name));
%% Generate cell array for visualization - scalar
scalar_stat_plot_cell = cell(3, 0);
% Global statistics
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.mask_volume_density, 'Vessel_volume_density', 'Volume density'};
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.mask_surface_area_density_mm2_mm3, 'Vessel_surface_area_density', 'Surface area density (mm^2/mm^3)'};
scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.link_length_density_m_mm3, 'Vessel_length_density', 'Length density (m/mm^3)'};

scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.capillary_volume_density, 'Capillary_volume_density', 'Volume fraction'};
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.capillary_length_density_m_mm3, 'Capillary_length_density', 'Length density (m/mm^3)'};
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.capillary_surface_area_density_mm2_mm3, 'Capillary_surface_area_density', 'Surface area density (mm^2/mm^3)'};

% capillary segment statistics
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.length, 'Median_capillary_length', 'Length(\mum)'};
% scalar_stat_plot_cell(:, end+1) = {cellfun(@(x) double(mean(x)), wb_pO2_stat_str.dt_lm.v), 'Average_local_DT_maximum', 'Distance (\mum)'};
% scalar_stat_plot_cell(:, end+1) = {wb_pO2_stat_str.local_dt_stat.mean, 'Average_tissue_vessel_distance', 'Distance (\mum)'};
% scalar_stat_plot_cell(:, end+1) = {cellfun(@(x) double(mean(x)), wb_pO2_stat_str.pO2_lm.v), 'Average_nromalized_tissue_pO2_local_minimum', 'Normalized pO_2'};

% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.cap2vsl_vol_fraction, 'Capillary-vessel_volume_ratio', 'Ratio'};
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.cap2vsl_length_fraction, 'Capillary-vessel_length_ratio', 'Ratio'};
% Loop
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.shortest_loop_geodesic_length, 'Average_number_of_edges_in_the_shortest_loop_of_capillary', 'Number of edges'};
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.shortest_loop_length, 'Average_length_of_the_shortest_loop_of_capillary', 'Length(\mum)'};
% Branching order
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.capillary_branching_order, 'Mean_capillary_branching_order', 'Capillary branching order'};
% Distance to large vessels
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.dist_to_nearest_noncapillary_mean, 'Medial_distance_to_nearest_noncapillary', 'Distance(\mum)'};
% Noncapillary
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.link_all_stat.mean.noncapillary_nearest_tissue_radius, 'Mean_noncapillary_nearest_tissue_radius' 'Nearest tissue radius(\mum)'};
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.link_all_stat.mean.noncapillary_nearest_capillary_num_vxl, 'Mean_nearest_caillary_total_vxl', 'Number of centerline voxel'};

% Node statistics
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.median.nearest_node_dist, 'Median_distance_to_the_nearest_node', 'Distance(\mum)'};
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.mean.path_to_nearest_neighbor_geodesic_length, 'Mean_geodesic_distance_to_nearest_node', 'Number of edges'};
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.num_data.degree ./ (0.240 * linear_scaling_factor) .^3, 'Node_density', 'Node density (mm^{-3})'};

% Segmentation quality
% scalar_stat_plot_cell(:, end+1) = {1 - wb_240_cube_stat_str.link_all_stat.mean.has_no_ep_Q, 'Fraction_of_vessel_with_endpoint', 'Fraction'};
% Radius
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.link_all_stat.median.dt_median, 'Vessel_median_radius', 'Radius (\mum)'};
% scalar_stat_plot_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.dt_median, 'Capillary_median_radius', 'Radius (\mum)'};
%% Generate cell array for visualization - anisotropy
anisotropy_plot_cell = cell(4, 0);
% Vessel 
anisotropy_plot_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_all_vw.fractional_anisotropy, ...
    'Volume_weighted_vessel_fractional_anisotropy', 'Fractional Anisotropy', wb_240_cube_stat_str.wb_ai_all_vw.svd_max_vec};
anisotropy_plot_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_all_vw.fa_z, ...
    'Volume_weighted_vessel_fractional_anisotropy_deviation_score', 'Fractional Anisotropy z-score', wb_240_cube_stat_str.wb_ai_all_vw.svd_max_vec};
% anisotropy_plot_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_all_vw.svd_value_ratio(:, 1), ...
%     'Volume_weighted_vessel_normalized_PC_value', 'Normalized principle value', wb_240_cube_stat_str.wb_ai_all_vw.svd_max_vec};

% Capillary 
anisotropy_plot_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_cap_vw.fractional_anisotropy, ...
    'Volume_weighted_capillary_fractional_anisotropy', 'Fractional Anisotropy', wb_240_cube_stat_str.wb_ai_cap_vw.svd_max_vec};
anisotropy_plot_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_cap_vw.fa_z, ...
    'Volume_weighted_capillary_fractional_anisotropy_deviation_score', 'Fractional Anisotropy z-score', wb_240_cube_stat_str.wb_ai_cap_vw.svd_max_vec};
% anisotropy_plot_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_cap_vw.svd_value_ratio(:, 1), ...
%     'Volume_weighted_capillary_normalized_PC_value', 'Normalized principle value', wb_240_cube_stat_str.wb_ai_cap_vw.svd_max_vec};
%% Generate videos to visualize the basic statistics on three planes
for vis_dir = 1 : 3
    proj_im_str = fun_vis_get_recon_max_proj_str(wb_max_proj_str, vis_dir);
    for iter_cell = 1 : size(scalar_stat_plot_cell, 2)
        if any(iter_cell == [3 5 7 8])
            save_layer_imageQ = true;
        else
            save_layer_imageQ = false;
        end        
        fun_vis_generate_local_stat_video(scalar_stat_plot_cell{:, iter_cell}, grid_info, proj_im_str, save_layer_imageQ, 5, 95);
    end
    
    for iter_cell = 1 : size(anisotropy_plot_cell, 2)
        save_layer_imageQ = true;
        fun_vis_generate_local_anisotropy_w_ori_video(anisotropy_plot_cell{:, iter_cell}, grid_info, proj_im_str, 'vec_proj', save_layer_imageQ, 5, 95);
    end
end