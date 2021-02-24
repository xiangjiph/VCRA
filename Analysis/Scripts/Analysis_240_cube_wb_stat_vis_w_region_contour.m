set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
% stack = 'ML20190124';
% stack = 'ML20200201';
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
anisotropy_data_name = 'Radius_dependence_anisotropy';
skel_version = '240_cube_rec';
vis_weight_method = 'volume_weighted';

recon_version = '240_cube_recon_sc';
plot_with_contour_Q = false;
%% Load data 
wb_max_proj_str = DataManager.load_analysis_data(dataset_name, stack, sprintf('%s_%s_%s_max_proj_data.mat', ...
    dataset_name, stack, mask_name), save_folder_name);
wb_240_cube_stat_str = DataManager.load_analysis_data(dataset_name, stack, sprintf('%s_%s_%s_240_cube_stat_data.mat', ...
    dataset_name, stack, mask_name), save_folder_name);
wb_pO2_stat_str = DataManager.load_analysis_data(dataset_name, stack, sprintf('%s_%s_%s_pO2_min_r_2um_sc_stat_data.mat', ...
    dataset_name, stack, skel_version));

rda_data = DataManager.load_data(DataManager.fp_analysis_data_file(...
    dataset_name, stack, sprintf('%s_%s_%s_%s_stat_data.mat', dataset_name, ...
    stack, skel_version, anisotropy_data_name)));
wb_vis_anisotropy_data = rda_data.(vis_weight_method);
%% Generate cell array for visualization - anisotropy
% Exclude cube that has less than 25% of the volume outside the brain 
is_vis_cube_Q = (wb_240_cube_stat_str.cube_in_brain_mask_ratio >= 0.1);
vis_radius_range_idx = 10;
vis_radius_min = rda_data.radius_min(vis_radius_range_idx);
vis_radius_max = rda_data.radius_max(vis_radius_range_idx);
vis_filename_template = sprintf('%s_%%s_r_%.2f_to_%.2f_um', vis_weight_method, ...
    vis_radius_min, vis_radius_max);
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
registration_vis_str.combined_regional_mask = combined_regional_mask;

region_name_table = struct;
region_name_table.Abbreviation = registration_vis_str.combined_region_abbrv';
region_name_table.Structure_name = vis_region_name';
region_name_table = struct2table(region_name_table);
tmp_fp = DataManager.fp_visualization_image(dataset_name, stack,  'video_region_name_LUT.csv', 'paper');
writetable(region_name_table, tmp_fp);
%% Generate cell array for visualization - scalar
plot_cell = cell(1, 0);
% Global statistics
% plot_cell{end + 1} = {wb_240_cube_stat_str.mask_volume_density, 'Vessel_volume_density', 'Volume density'};
% plot_cell{end + 1} = {wb_240_cube_stat_str.mask_surface_area_density_mm2_mm3, 'Vessel_surface_area_density', 'Surface area density (mm^2/mm^3)'};
plot_cell{end+1} = {wb_240_cube_stat_str.link_length_density_m_mm3, 'Vessel_length_density', 'Vessel length density (m/mm^3)'};

% plot_cell{end + 1} = {wb_240_cube_stat_str.capillary_volume_density, 'Capillary_volume_density', 'Volume fraction'};
plot_cell{end+1} = {wb_240_cube_stat_str.capillary_length_density_m_mm3, 'Capillary_length_density', 'Capillary length density (m/mm^3)'};
% plot_cell{end + 1} = {wb_240_cube_stat_str.capillary_surface_area_density_mm2_mm3, 'Capillary_surface_area_density', 'Surface area density (mm^2/mm^3)'};

% capillary segment statistics
% plot_cell{end + 1} = {wb_240_cube_stat_str.link_cap_stat.median.length, 'Median_capillary_length', 'Length(\mum)'};
% plot_cell{end + 1} = {cellfun(@(x) double(mean(x)), wb_pO2_stat_str.dt_lm.v), 'Average_local_DT_maximum', 'Distance (\mum)'};
% % plot_cell{end + 1} = {wb_pO2_stat_str.local_dt_stat.mean, 'Average_tissue_vessel_distance', 'Distance (\mum)'};
% plot_cell{end + 1} = {cellfun(@(x) double(mean(x)), wb_pO2_stat_str.pO2_lm.v), 'Average_nromalized_tissue_pO2_local_minimum', 'Normalized pO_2'};
% Loop
% plot_cell{end + 1} = {wb_240_cube_stat_str.link_cap_stat.mean.shortest_loop_geodesic_length, 'Average_number_of_edges_in_the_shortest_loop_of_capillary', 'Number of edges'};
% plot_cell{end + 1} = {wb_240_cube_stat_str.link_cap_stat.mean.shortest_loop_length, 'Average_length_of_the_shortest_loop_of_capillary', 'Length(\mum)'};
% Branching order
% plot_cell{end + 1} = {wb_240_cube_stat_str.link_cap_stat.mean.capillary_branching_order, 'Mean_capillary_branching_order', 'Capillary branching order'};
% Distance to large vessels
% plot_cell{end + 1} = {wb_240_cube_stat_str.link_cap_stat.median.dist_to_nearest_noncapillary_mean, 'Medial_distance_to_nearest_noncapillary', 'Distance(\mum)'};
% Noncapillary
% plot_cell{end + 1} = {wb_240_cube_stat_str.link_all_stat.mean.noncapillary_nearest_tissue_radius, 'Mean_noncapillary_nearest_tissue_radius' 'Nearest tissue radius(\mum)'};
% plot_cell{end + 1} = {wb_240_cube_stat_str.link_all_stat.mean.noncapillary_nearest_capillary_num_vxl, 'Mean_nearest_caillary_total_vxl', 'Number of centerline voxel'};

% Node statistics
% plot_cell{end + 1} = {wb_240_cube_stat_str.node_stat.median.nearest_node_dist, 'Median_distance_to_the_nearest_node', 'Distance(\mum)'};
% plot_cell{end + 1} = {wb_240_cube_stat_str.node_stat.mean.path_to_nearest_neighbor_geodesic_length, 'Mean_geodesic_distance_to_nearest_node', 'Number of edges'};
% plot_cell{end + 1} = {wb_240_cube_stat_str.node_stat.num_data.degree ./ (0.240 * linear_scaling_factor) .^3, 'Node_density', 'Node density (mm^{-3})'};

% Segmentation quality
% plot_cell{end + 1} = {1 - wb_240_cube_stat_str.link_all_stat.mean.has_no_ep_Q, 'Fraction_of_vessel_with_endpoint', 'Fraction'};
%% Generate cell array for visualization - anisotropy
% plot_cell{end + 1} = {wb_240_cube_stat_str.wb_ai_cap_lw.fa_z, ...
%     'Length_weighted_capillary_fractional_anisotropy_deviation_score', 'Fractional Anisotropy z-score', wb_240_cube_stat_str.wb_ai_cap_lw.svd_max_vec};
% plot_cell{end + 1} = {wb_240_cube_stat_str.wb_ai_cap_lw.fractional_anisotropy, ...
%     'Length_weighted_capillary_fractional_anisotropy', 'Fractional Anisotropy', wb_240_cube_stat_str.wb_ai_cap_lw.svd_max_vec};
% plot_cell{end + 1} = {wb_240_cube_stat_str.wb_ai_cap_vw.svd_value_ratio(:, 1), ...
%     'Volume_weighted_capillary_normalized_PC_value', 'Normalized principle value', wb_240_cube_stat_str.wb_ai_cap_vw.svd_max_vec};

plot_cell{end + 1} = {wb_vis_anisotropy_data.fa(:, vis_radius_range_idx), ...
    sprintf(vis_filename_template, 'FA'), 'Fractional Anisotropy', wb_vis_anisotropy_data.svd_max_vec(:, :, vis_radius_range_idx)};

plot_cell{end + 1} = {wb_vis_anisotropy_data.fa_z(:, vis_radius_range_idx), ...
    sprintf(vis_filename_template, 'FAz'), 'FA_z', wb_vis_anisotropy_data.svd_max_vec(:, :, vis_radius_range_idx)};

fa_p = wb_vis_anisotropy_data.fa_p(:, vis_radius_range_idx);
fa_p(fa_p == 0) = 5e-5;
fa_p_log = -log10(fa_p);
plot_cell{end + 1} = {fa_p_log, ...
    sprintf(vis_filename_template, 'FAp'), '- log_{10}(FA_p)', wb_vis_anisotropy_data.svd_max_vec(:, :, vis_radius_range_idx)};
%% Generate videos to visualize the basic statistics on three planes
for vis_dir = 1 : 3
    proj_im_str = fun_vis_get_recon_max_proj_str(wb_max_proj_str, vis_dir);
    
    if plot_with_contour_Q
        tmp_registration_vis_str = fun_vis_permute_registration_str(registration_vis_str, vis_dir);
    else
        tmp_registration_vis_str = [];
    end
    
    for iter_cell = 1 : numel(plot_cell)
%         if any(iter_cell == [1 2 3 4])
%             save_layer_imageQ = true;
%         else
            save_layer_imageQ = false;
%         end
        tmp_plot_data_cell = plot_cell{iter_cell};
%         tmp_plot_data_cell{1}(~is_vis_cube_Q) = nan;
        switch numel(tmp_plot_data_cell)
            case 3
                fun_vis_generate_local_stat_video_w_region_boundary(tmp_plot_data_cell{:}, ...
                    grid_info, proj_im_str, tmp_registration_vis_str, save_layer_imageQ, linear_scaling_factor, 5, 95);
            case 4
                fun_vis_generate_local_anisotropy_w_ori_video_w_boundary(tmp_plot_data_cell{:}, ...
                    grid_info, proj_im_str, tmp_registration_vis_str, ...
                    'vec_proj', save_layer_imageQ, linear_scaling_factor, 5, 95);
        end
    end
end