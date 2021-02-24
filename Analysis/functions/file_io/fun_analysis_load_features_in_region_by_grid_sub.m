function [features, varargout] = fun_analysis_load_features_in_region_by_grid_sub(dataset_name, stack, mask_version, grid_sub)

persistent DataManager;
if isempty(DataManager)
    DataManager = FileManager;
end

if nargout > 1
    compute_region_stat_Q = true;
else
    compute_region_stat_Q = false;
end

num_cube = size(grid_sub, 1);
report_cube_ind = round(num_cube / 10);
link_feature = cell(num_cube, 1);
node_feature = cell(num_cube, 1);
tmp_tic = tic;

if compute_region_stat_Q
    fprintf('Compute the 240-cube statistics\n');
    subgrid_stat = struct;
    subgrid_stat.grid_sub = grid_sub';
    
    [subgrid_stat.link_all, subgrid_stat.link_cap, subgrid_stat.node]  = deal(cell(num_cube, 1));
    
    [subgrid_stat.surface_area_mm2mm3, subgrid_stat.volume_ratio, ...
        subgrid_stat.length_density_m_mm3, subgrid_stat.length_density_cap_m_mm3, ...
        subgrid_stat.surface_area_cap_mm2mm3, subgrid_stat.volume_ratio_cap] = deal(nan(num_cube, 1));
    
    [subgrid_stat.ai_all_vw_fa, subgrid_stat.ai_all_vw_fa_z, subgrid_stat.ai_all_vw_min2max_z, subgrid_stat.ai_all_vw_svd1_z, ...
        subgrid_stat.ai_cap_vw_fa, subgrid_stat.ai_cap_vw_fa_z, subgrid_stat.ai_cap_vw_min2max_z, subgrid_stat.ai_cap_vw_svd1_z, ...
        subgrid_stat.ai_all_vw_svrmax, subgrid_stat.ai_cap_vw_svrmax] = deal(nan(num_cube, 1));
    
    [subgrid_stat.ai_all_vw_vec, subgrid_stat.ai_cap_vw_vec, subgrid_stat.ai_noncap_vw_vec] = deal(nan(3, num_cube));
    
    capillary_max_r = 3.5;
end

for iter_bbox = 1 : num_cube
   tmp_idx_1 = grid_sub(iter_bbox, 1);
   tmp_idx_2 = grid_sub(iter_bbox, 2);
   tmp_idx_3 = grid_sub(iter_bbox, 3);
   try 
      tmp_mask_str = DataManager.load_block_mask(dataset_name, stack, mask_version, ...
          tmp_idx_1, tmp_idx_2, tmp_idx_3);
      tmp_mask_str = fun_analysis_delete_unused_features_in_recon_mask_str(tmp_mask_str);
      link_feature{iter_bbox} = tmp_mask_str.link.features;
      node_feature{iter_bbox} = tmp_mask_str.node.features;       
      
      if compute_region_stat_Q          
          % Compute basic feature statistics
          subgrid_stat.link_all{iter_bbox} = fun_analysis_get_graph_feature_basic_stat(tmp_mask_str.link.features);
          
          subgrid_stat.link_cap{iter_bbox} = fun_analysis_get_graph_feature_basic_stat(tmp_mask_str.link.features,...
              tmp_mask_str.link.features.dt_median <= capillary_max_r);
          
          subgrid_stat.node{iter_bbox} = fun_analysis_get_graph_feature_basic_stat(tmp_mask_str.node.features);
          % Compute the local orientation of the noncapillaries
          if any(tmp_mask_str.link.features.is_large_vessel_Q)
              tmp_non_cap_ep2ep_vec_vw = bsxfun(@times, ...
                  tmp_mask_str.link.features.ep1_to_ep2_direction_vec(tmp_mask_str.link.features.is_large_vessel_Q, :), ...
                  tmp_mask_str.link.features.in_bbox_volume(tmp_mask_str.link.features.is_large_vessel_Q));
              tmp_is_nan_Q = any(isnan(tmp_non_cap_ep2ep_vec_vw), 2);
              tmp_non_cap_ep2ep_vec_vw = tmp_non_cap_ep2ep_vec_vw(~tmp_is_nan_Q, :);
              if ~isempty(tmp_non_cap_ep2ep_vec_vw)
                  [svd_U, ~, ~] = svd(cov(cat(1, tmp_non_cap_ep2ep_vec_vw, -tmp_non_cap_ep2ep_vec_vw)), 'econ');
                  subgrid_stat.ai_noncap_vw_vec(:, iter_bbox) = svd_U(:, 1);
              end
          end
          
          cube_volume_mm3 = prod(tmp_mask_str.block_size / 1e3);
          subgrid_stat.surface_area_mm2mm3(iter_bbox) = tmp_mask_str.mask_surface_area / (1e6 * cube_volume_mm3);
          subgrid_stat.volume_ratio(iter_bbox) = tmp_mask_str.mask_volume_ratio;          
          subgrid_stat.length_density_m_mm3(iter_bbox) = tmp_mask_str.total_length_um / (1e6 * cube_volume_mm3);
          subgrid_stat.length_density_cap_m_mm3(iter_bbox) = tmp_mask_str.total_capillary_length_um / (1e6 * cube_volume_mm3);
          subgrid_stat.surface_area_cap_mm2mm3(iter_bbox) = tmp_mask_str.total_capillary_surface_area_um2 / (1e6 * cube_volume_mm3);
          subgrid_stat.volume_ratio_cap(iter_bbox) = tmp_mask_str.total_capillary_volume_um3 / (1e9 * cube_volume_mm3);

          if isfield(tmp_mask_str.anisotropy_all_vw, 'fractional_anisotropy') &&...
                  ~isempty(tmp_mask_str.anisotropy_all_vw.fractional_anisotropy)
              subgrid_stat.ai_all_vw_fa(iter_bbox) = tmp_mask_str.anisotropy_all_vw.fractional_anisotropy;
              subgrid_stat.ai_all_vw_fa_z(iter_bbox) = tmp_mask_str.anisotropy_all_vw.fa_z;
              
              subgrid_stat.ai_all_vw_min2max_z(iter_bbox) = tmp_mask_str.anisotropy_all_vw.min2max_z;
              subgrid_stat.ai_all_vw_svd1_z(iter_bbox) = tmp_mask_str.anisotropy_all_vw.svd_1_z;
              
              subgrid_stat.ai_all_vw_vec(:, iter_bbox) = tmp_mask_str.anisotropy_all_vw.svd_max_vec;
              subgrid_stat.ai_all_vw_svrmax(iter_bbox) = tmp_mask_str.anisotropy_all_vw.svd_value_ratio(1);
          end
          
          if isfield(tmp_mask_str.anisotropy_capillary_vw, 'fractional_anisotropy') &&...
                  ~isempty(tmp_mask_str.anisotropy_capillary_vw.fractional_anisotropy)
              subgrid_stat.ai_cap_vw_fa(iter_bbox) = tmp_mask_str.anisotropy_capillary_vw.fractional_anisotropy;
              subgrid_stat.ai_cap_vw_fa_z(iter_bbox) = tmp_mask_str.anisotropy_capillary_vw.fa_z;
              
              subgrid_stat.ai_cap_vw_min2max_z(iter_bbox) = tmp_mask_str.anisotropy_capillary_vw.min2max_z;
              subgrid_stat.ai_cap_vw_svd1_z(iter_bbox) = tmp_mask_str.anisotropy_capillary_vw.svd_1_z;
              subgrid_stat.ai_cap_vw_vec(:, iter_bbox) = tmp_mask_str.anisotropy_capillary_vw.svd_max_vec;              
              subgrid_stat.ai_cap_vw_svrmax(iter_bbox) = tmp_mask_str.anisotropy_capillary_vw.svd_value_ratio(1);
          end
      end
   catch ME
       fprintf('Failed to load the mask of block (%d, %d, %d)\n', tmp_idx_1, tmp_idx_2, tmp_idx_3);
       fun_rethrow_error_message_without_error(ME);
   end
   
   if mod(iter_bbox, report_cube_ind) == 0
       fprintf('Finish loading %d / %d. Elapsed time is %f seconds\n', iter_bbox, num_cube, toc(tmp_tic));
   end   
end
% Remove empty tables
link_feature = link_feature(cellfun(@(x) size(x, 1), link_feature) > 0);
node_feature = node_feature(cellfun(@(x) size(x, 1), node_feature) > 0);
features.link = cat(1, link_feature{:});
features.node = cat(1, node_feature{:});
assert(istable(features.link) && istable(features.node));

if compute_region_stat_Q
    varargout{1} = subgrid_stat;
end

end
%% Subfunction
function recon_mask_str = fun_analysis_delete_unused_features_in_recon_mask_str(recon_mask_str)

link_feature_to_rm = {'link_com', 'ep1_sub', 'ep2_sub', 'mid_sub', 'dt_ep1', 'dt_ep2', ...
    'dt_min', 'dt_max', 'dt_std', 'ep2ep_angle_azimuth_deg', 'ep2ep_angle_elevation_deg', ...
    'ep1_direction_vec', 'ep2_direction_vec', 'cc_sub_pca1_vec', 'cc_sub_pca2_vec', ...
    'cc_sub_pca3_vec', 'cc_sub_cov_eig_val', 'shortest_loop_link_label', 'shortest_loop_node_label'};

node_feature_to_rm = {'nearest_node_label', 'center_of_mass', 'nearest_tissue_dt_max', ...
    'nearest_tissue_dt_median', 'nearest_tissue_dt_mean', 'nearest_tissue_volume', 'neighbor_node_label'};
assert(istable(recon_mask_str.link.features))
recon_mask_str.link.features = removevars(recon_mask_str.link.features, ...
    link_feature_to_rm);
assert(istable(recon_mask_str.node.features))
recon_mask_str.node.features = removevars(recon_mask_str.node.features, ...
    node_feature_to_rm);
end
