function recon_str = fun_analysis_sc_reconstruction_result(recon_str, linear_scaling_factor, radius_scaling_factor)


% Vessel radius is not scaled at the moment. 
%% Scale link features
link_linear_scale_features = {'length', 'ep2ep_dist', 'shortest_loop_length', ...
    'nearest_tissue_dt_mean', 'nearest_tissue_dt_max', 'nearest_tissue_dt_median', ...
    'dist_rm_lk_2_nlm', 'lk_dt_af_rm_mean', ...
    'lk_dt_af_rm_median', 'lk_dt_af_rm_max', 'pt_vol_dt_max', 'dist_rm_lk_2_pt_max_dt', ...
    'pt_vol_dt_mean', 'pt_vol_dt_median', 'nearest_tissue_radius', 'diff_dp_d_mean'...
    'noncapillary_nearest_tissue_dt_mean', 'noncapillary_nearest_tissue_dt_median', ...
    'noncapillary_nearest_tissue_dt_max', 'noncapillary_nearest_capillary_num_vxl', ...
    'dist_ep1_to_nearest_noncapillary', 'dist_ep2_to_nearest_noncapillary', 'dist_to_nearest_noncapillary_mean', ...
    'dist_to_nearest_noncapillary_median', 'dist_to_brain_surface_mean', 'in_bbox_length', ...
    'in_bbox_surface_area', 'surface_area', 'in_bbox_volume', 'volume'};

link_cubic_scale_features = {'nearest_tissue_volume'};

link_r_scale_features = {'in_bbox_surface_area', 'surface_area'};
link_r_sqred_scale_features = {'in_bbox_volume', 'volume'};

assert(istable(recon_str.link.features));
if istable(recon_str.link.features)
    link_feature_name = recon_str.link.features.Properties.VariableNames;
elseif isstruct(recon_str.link.features)
    link_feature_name = fieldnames(recon_str.link.features);
else
    error('Unrecognized feature data');
end

recon_str.link.features = fun_analysis_sc_structure_fields(recon_str.link.features, ...
    intersect(link_linear_scale_features, link_feature_name), linear_scaling_factor);

recon_str.link.features = fun_analysis_sc_structure_fields(recon_str.link.features, ...
    intersect(link_cubic_scale_features, link_feature_name), linear_scaling_factor ^ 3);

recon_str.link.features = fun_analysis_sc_structure_fields(recon_str.link.features, ...
    intersect(link_r_scale_features, link_feature_name), radius_scaling_factor);
% 
recon_str.link.features = fun_analysis_sc_structure_fields(recon_str.link.features, ...
    intersect(link_r_sqred_scale_features, link_feature_name), radius_scaling_factor ^ 2);
%% Scale node features
node_linear_scale_features = {'nearest_node_dist', 'radial_node_count_edge', 'neighbor_node_dist'...
    'link_length', 'link_length_max', 'link_length_min', 'link_length_mean', 'link_length_median', ...
    'path_to_neighbor_length', 'path_to_nearest_neighbor_length', 'dist_to_brain_surface_mean'};
node_neg_cube_scale_features = {'radial_node_density'};
assert(istable(recon_str.node.features));
if istable(recon_str.node.features)
    node_feature_name = recon_str.node.features.Properties.VariableNames;
elseif isstruct(recon_str.link.features)
    node_feature_name = fieldnames(recon_str.node.features);
else
    error('Unrecognized feature data');
end

recon_str.node.features = fun_analysis_sc_structure_fields(recon_str.node.features, ...
    intersect(node_linear_scale_features, node_feature_name), linear_scaling_factor);

recon_str.node.features = fun_analysis_sc_structure_fields(recon_str.node.features, ...
    intersect(node_neg_cube_scale_features, node_feature_name), linear_scaling_factor ^ (-3));
%% Scale statistics
% The logic is to assume that vessel is a cylinder and therefore the length and
% the radius can be scaled separately
recon_str.stat = fun_analysis_sc_structure_fields(recon_str.stat, ...
    {'mask_surface_area_density_mm2_mm3', 'mask_volume_density', ...
    'link_surface_area_density_mm2_mm3', 'link_volume_density', 'link_length_density_m_mm3', ...
    'capillary_surface_area_density_mm2_mm3', 'capillary_volume_density', 'capillary_length_density_m_mm3'}, linear_scaling_factor^(-2));

recon_str.stat = fun_analysis_sc_structure_fields(recon_str.stat, ...
    {'mask_surface_area_density_mm2_mm3', 'link_surface_area_density_mm2_mm3', 'capillary_surface_area_density_mm2_mm3'}, radius_scaling_factor);

recon_str.stat = fun_analysis_sc_structure_fields(recon_str.stat, ...
    {'capillary_volume_density', 'link_volume_density', 'mask_volume_density'}, radius_scaling_factor ^ 2);
end