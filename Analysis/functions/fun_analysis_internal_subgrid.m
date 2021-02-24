function exit_code = fun_analysis_internal_subgrid(dataset_name, stack, load_skl_name, grid_c_label, opt)
% fun_analysis_internal_subgrid is a wrap up function modified from
% Analysis_generate_combined_5_o_2_graph. This function load the combined
% grid skeleton, generate the graph, analyze the graph features and write
% the internal subgrid statistics to the hard drive

% To do list: 
% 1. Add distance transform of the annotated mask. Compute the distance
% between the link and node cc to the surface of the brain. The distance
% transform should be computed once and made persistent for the same task. 
%% Parameters
if nargin < 5
    % This part need to be replaced by loading the default setting from the
    % XML file. 
    DataManager = FileManager;
%     grid_c_info = DataManager.load_grid(dataset_name, stack, grid_c_version);
%     opt_comp_graph_feature = struct;
%     computeQ = struct;
%     computeQ.basic_link_feature = true;
%     computeQ.basic_node_feature = true;
%     computeQ.basic_reconstruction = true;
%     computeQ.node_path_to_nearest_neighbor = true;
%     computeQ.link_shortest_path = true;
%     computeQ.dist_tissue_2_vessel = true;
%     computeQ.dimension = false;
%     computeQ.max_z_proj = true;
%     opt_comp_graph_feature.computeQ = computeQ;
%     opt_comp_graph_feature.vis_dim_fit = false;
%     opt_comp_graph_feature.dim_fit_cutoff_length = 25;
%     opt_comp_graph_feature.recon_max_error_rate = 0.1;
%     opt_comp_graph_feature.merge_neighbor_nodes_Q = false;
%     % Option for local statistics
%     opt_comp_graph_feature.capillary_max_radius = 3.5; % um;
%     overwrite_Q = false;
%     subgrid_mask_name = sprintf('%s_recon', skl_grid_name);
else
    if ~isa(opt, 'struct')
        if isfile(opt)
            opt = load(opt);
        else
            error('The input opt should be either a valid path to the mat file or a MATLAB structure');
        end
    end
    DataManager = opt.DataManager;
    grid_c_info = opt.grid_c_info;
    opt_comp_graph_feature = opt.graph_feature;
    overwrite_Q = opt.overwrite_Q;
    output_mask_name = opt.output_mask_name;
    output_graph_name = opt.output_graph_name;
end
grid_info = grid_c_info.grid_ori;
%%
tmp_grid_c_sub = grid_c_info.bbox_grid_sub_list(grid_c_label, :);
tmp_grid_c_ind = sub2ind(grid_c_info.grid_size, tmp_grid_c_sub(1), ...
    tmp_grid_c_sub(2), tmp_grid_c_sub(3));
assert(grid_c_info.bbox_xyz_label_array(tmp_grid_c_ind) == grid_c_label, 'Incorrect grid ind');
% Number of valid subgrid in the grid_c
tmp_subgrid_valid_idx = grid_c_info.internal_subgrid_valid_idx{tmp_grid_c_ind};
if isempty(tmp_subgrid_valid_idx)
    fprintf('No valid internal sub-grid in this combined grid\n');
    exit_code = 2;
    return
end
%% Combined grid information
% graph_fp = DataManager.fp_graph_in_block_file(dataset_name, stack, output_graph_name, ...
%     tmp_grid_c_sub(1), tmp_grid_c_sub(2), tmp_grid_c_sub(3));
% if isfile(graph_fp) && ~overwrite_Q
%     fprintf('Graph alreadly exist. Do not overwrite\n');
%     exit_code = 3;
%     return;
% end
data_info = fun_graph_get_grid_c_block_info(grid_c_info, grid_c_label);
%% Load all the subgrid data
grid_c_subgrid_sub = grid_c_info.sub_grid_sub{tmp_grid_c_ind};
grid_c_subgrid_sub_min = grid_c_info.bbox_xyz_mmll_grid_list(grid_c_label, 1:3);
grid_c_subgrid_sub_local = bsxfun(@minus, grid_c_subgrid_sub, grid_c_subgrid_sub_min - 1);

grid_c_subgrid_valid_array = grid_c_info.sub_grid_label_array{tmp_grid_c_ind} > 0;

grid_c_num_subgrid = size(grid_c_subgrid_sub, 1);
assert(nnz(grid_c_subgrid_valid_array) == grid_c_num_subgrid, 'Number of valid subgrid does not match');
grid_c_mmxx_pixel = grid_c_info.bbox_xyz_mmxx_pixel_list(grid_c_label, :);
grid_c_range_pixel = grid_c_mmxx_pixel(4:6) - grid_c_mmxx_pixel(1:3) + 1;
% Load skeleton and combined them
grid_c_skl_ind = cell(grid_c_num_subgrid, 1);
grid_c_r = cell(grid_c_num_subgrid, 1);
grid_c_bbox = zeros(6, grid_c_num_subgrid);
for iter_cube = 1 : grid_c_num_subgrid
    tmp_idx_1 = grid_c_subgrid_sub(iter_cube, 1);
    tmp_idx_2 = grid_c_subgrid_sub(iter_cube, 2);
    tmp_layer = grid_c_subgrid_sub(iter_cube, 3);
    try
        tmp_skel = DataManager.load_block_skl(dataset_name, stack, load_skl_name, tmp_idx_1, tmp_idx_2, tmp_layer);
    catch
        fprintf('Fail to read the skeleton file of cube (%d, %d, %d) . Skip this block\n', tmp_idx_1, tmp_idx_2, tmp_layer);
        continue;
    end
    tmp_bbox_mmxx = fun_get_block_valid_bbox_mmxx(grid_c_subgrid_sub_local(iter_cube, :),...
        grid_c_subgrid_valid_array, tmp_skel.block_size, grid_info.block_overlap);
    tmp_skl_sub = fun_ind2sub(tmp_skel.block_size, tmp_skel.ind);
    tmp_in_bbox_Q = all( bsxfun(@ge, tmp_skl_sub, tmp_bbox_mmxx(1:3)), 2) & ...
        all( bsxfun(@le, tmp_skl_sub, tmp_bbox_mmxx(4:6)), 2);
    tmp_skl_sub_valid = tmp_skl_sub(tmp_in_bbox_Q , : );
    % Convert to coordinate in the combined grid
    tmp_cube_disp = tmp_skel.global_bbox_mmll(1:3) - grid_c_mmxx_pixel(1:3);
    tmp_skl_sub_in_grid_c = tmp_skl_sub_valid + tmp_cube_disp;
    
    grid_c_bbox(:, iter_cube) = tmp_bbox_mmxx + [tmp_cube_disp, tmp_cube_disp];
    
    grid_c_skl_ind{iter_cube} = sub2ind(grid_c_range_pixel, tmp_skl_sub_in_grid_c(:,1), ...
        tmp_skl_sub_in_grid_c(:, 2), tmp_skl_sub_in_grid_c(:,3));
    grid_c_r{iter_cube} = tmp_skel.r(tmp_in_bbox_Q);
end
grid_c_skl_ind = cat(1, grid_c_skl_ind{:});
grid_c_r = cat(1, grid_c_r{:});
assert(numel(grid_c_skl_ind) == numel(grid_c_r), 'Number of skeleton voxels does not equal to the number of radius');
% For debug
tmp_num_skl_ind = numel(grid_c_skl_ind);
% It's possible to have duplicated skeleton voxels. One case is when 6
% neighbors are all empty but the remaining 20 neighbors of the 26
% neighbors are not empty.
[grid_c_skl_ind, tmp_unique_idx, ~] = unique(grid_c_skl_ind);
grid_c_r = grid_c_r(tmp_unique_idx);

tmp_num_skl_ind_unique = numel(grid_c_skl_ind);
if tmp_num_skl_ind ~= tmp_num_skl_ind_unique
    if all(grid_c_subgrid_valid_array(:))
        error('Duplicated skeleton voxel exist. Current block index is %d\n', grid_c_label);
    else
        fprintf('Exist duplicated skeleton voxels. Duplicated ratio is %f\n', (tmp_num_skl_ind - tmp_num_skl_ind_unique) / tmp_num_skl_ind_unique);
    end
end
if tmp_num_skl_ind_unique == 0
   fprintf('Empty block! Return\n');
   exit_code = 4;
   return;
end
% Convert to graph
vessel_graph = fun_skeleton_to_graph(grid_c_skl_ind, grid_c_range_pixel);
vessel_graph.radius = sparse(grid_c_skl_ind, ones(vessel_graph.num.skeleton_voxel, 1), ...
    double(grid_c_r), vessel_graph.num.block_voxel, 1);
% Merge neighbor nodes
% if isfield(opt_graph_feature, 'merge_neighbor_nodes_Q') && opt_graph_feature.merge_neighbor_nodes_Q
%     link_label_to_merge = fun_graph_get_link_label_to_be_merged(vessel_graph.link.cc_ind,...
%         vessel_graph.radius, opt_graph_feature.max_merge_length, opt_graph_feature.always_merge_max_length);
%     vessel_graph = fun_graph_merge_node_by_link_label(vessel_graph, merged_link_label, true);
% end
%% Get the links and nodes that have at least one voxel in the sub-grid bounding box
tmp_internal_subgrid_bbox_mmxx_dataset = grid_c_info.internal_subgrid_bbox_mmxx{tmp_grid_c_ind};
tmp_internal_subgrid_bbox_mmxx_grid_c = tmp_internal_subgrid_bbox_mmxx_dataset - ...
    [grid_c_mmxx_pixel(1:3) - 1, grid_c_mmxx_pixel(1:3) - 1];
assert(~iscolumn(tmp_internal_subgrid_bbox_mmxx_grid_c), 'The internal subgrid bbox mmxx list is a column vector');
tmp_num_internal_subgrid = size(tmp_internal_subgrid_bbox_mmxx_grid_c, 1);

tmp_link_label_in_bbox_cell = cell(tmp_num_internal_subgrid, 1);
tmp_node_label_in_bbox_cell = cell(tmp_num_internal_subgrid, 1);
tmp_link_voxel_in_bbox_ratio_cell = cell(tmp_num_internal_subgrid, 1);
for iter_int_grid = 1 : tmp_num_internal_subgrid
    tmp_int_bbox = tmp_internal_subgrid_bbox_mmxx_grid_c(iter_int_grid, :);
    % All the links that pass through the bounding box
    tmp_link_sub = fun_ind2sub(vessel_graph.num.mask_size, vessel_graph.link.pos_ind);
    tmp_node_sub = fun_ind2sub(vessel_graph.num.mask_size, vessel_graph.node.pos_ind);
    tmp_link_voxel_in_bbox_Q = all(bsxfun(@ge, tmp_link_sub, tmp_int_bbox(1:3)), 2) & ...
        all(bsxfun(@le, tmp_link_sub, tmp_int_bbox(4:6)), 2);
    tmp_node_voxel_in_bbox_Q = all(bsxfun(@ge, tmp_node_sub, tmp_int_bbox(1:3)), 2) & ...
        all(bsxfun(@le, tmp_node_sub, tmp_int_bbox(4:6)), 2);
    tmp_voxel_link_label = full(vessel_graph.link.map_ind_2_label(...
        vessel_graph.link.pos_ind(tmp_link_voxel_in_bbox_Q)));
    tmp_voxel_node_label = full(vessel_graph.node.map_ind_2_label(...
        vessel_graph.node.pos_ind(tmp_node_voxel_in_bbox_Q)));
    
    [tmp_voxel_link_label_cell, tmp_link_in_bbox_label] = fun_bin_data_to_idx_list(tmp_voxel_link_label);
    
    tmp_node_in_bbox_label = unique(tmp_voxel_node_label);
    
    tmp_link_num_voxel = vessel_graph.link.num_voxel_per_cc(tmp_link_in_bbox_label);
    tmp_link_num_voxel_in_bbox = cellfun(@numel, tmp_voxel_link_label_cell);
    tmp_link_in_bbox_ratio = (tmp_link_num_voxel_in_bbox ./ tmp_link_num_voxel);
    
    tmp_link_label_in_bbox_cell{iter_int_grid} = tmp_link_in_bbox_label;
    tmp_node_label_in_bbox_cell{iter_int_grid} = tmp_node_in_bbox_label;
    tmp_link_voxel_in_bbox_ratio_cell{iter_int_grid} = tmp_link_in_bbox_ratio;
end
tmp_internal_link_lable = unique(cat(1, tmp_link_label_in_bbox_cell{:}));
tmp_num_internal_link = numel(tmp_internal_link_lable);
tmp_internal_node_label = unique(cat(1, tmp_node_label_in_bbox_cell{:}));
tmp_num_internal_node = numel(tmp_internal_node_label);
tmp_link_label_2_list_ind = sparse(tmp_internal_link_lable, ones(tmp_num_internal_link, 1), ...
    1 : tmp_num_internal_link);
tmp_node_label_2_list_ind = sparse(tmp_internal_node_label, ones(tmp_num_internal_node, 1), ...
    1 : tmp_num_internal_node);
%% Compute graph features for internal node and links
%     disp('Compute internal links and nodes features');]
% Handle the error message outside the function
% try
if isempty(tmp_internal_link_lable) && isempty(tmp_internal_node_label)
    fprintf('No links in the internal subgrids\n');
    exit_code = 5;
    return;
end
vessel_graph.info = data_info;
[vessel_graph, recon_mask] = fun_analysis_compute_graph_features_by_labels(vessel_graph, tmp_internal_link_lable, tmp_internal_node_label, opt_comp_graph_feature);

if ~isempty(tmp_internal_link_lable)
    assert(size(vessel_graph.link.features, 1) == numel(tmp_internal_link_lable), 'Number of items in the table does not equal the number of link labels');
    vessel_graph.link.features.has_1_ep_Q = ~all(vessel_graph.link.connected_node_label(tmp_internal_link_lable, :), 2) & ...
        any(vessel_graph.link.connected_node_label(tmp_internal_link_lable, :), 2);
    vessel_graph.link.features.has_no_ep_Q = all(vessel_graph.link.connected_node_label(tmp_internal_link_lable, :), 2);
end

%% Output
% Divide the computed features into internal 240 cubes and save the result
tmp_internal_subgrid_bbox_mmll_dataset = grid_c_info.internal_subgrid_bbox_mmll{tmp_grid_c_ind};
tmp_internal_subgrid_bbox_mmll_grid_c = tmp_internal_subgrid_bbox_mmll_dataset;
tmp_internal_subgrid_bbox_mmll_grid_c(:, 1:3) = tmp_internal_subgrid_bbox_mmll_grid_c(:, 1:3) - (grid_c_mmxx_pixel(1:3) - 1);
tmp_internal_subgrid_sub = grid_c_info.internal_subgrid_sub{tmp_grid_c_ind};

grid_c_skl_sub = fun_ind2sub(vessel_graph.num.mask_size, grid_c_skl_ind);
grid_c_skl_global_sub = grid_c_skl_sub + grid_c_mmxx_pixel(1:3) - 1;
for iter_int_grid = 1 : tmp_num_internal_subgrid
    tmp_subgrid_graph = [];
    tmp_subgrid_recon = struct;
    % Add subgrid information
    tmp_subgrid_recon.dataset_name = dataset_name;
    tmp_subgrid_recon.stack = stack;
    tmp_subgrid_recon.version = output_mask_name;
    tmp_subgrid_recon.grid_version = grid_info.version;
    tmp_subgrid_recon.grid_idx_1 = tmp_internal_subgrid_sub(iter_int_grid, 1);
    tmp_subgrid_recon.grid_idx_2 = tmp_internal_subgrid_sub(iter_int_grid, 2);
    tmp_subgrid_recon.grid_layer = tmp_internal_subgrid_sub(iter_int_grid, 3);
    tmp_subgrid_recon.global_bbox_mmxx = tmp_internal_subgrid_bbox_mmxx_dataset(iter_int_grid, :);
    tmp_subgrid_recon.global_bbox_mmll = tmp_internal_subgrid_bbox_mmll_dataset(iter_int_grid, :);
    tmp_subgrid_recon.global_block_size = grid_info.data_size;
    tmp_subgrid_recon.block_size = tmp_subgrid_recon.global_bbox_mmll(4:6);
    % Graph 
    tmp_skl_vxl_in_subgrid_Q = fun_voxel_sub_in_bbox_mmxx_Q(grid_c_skl_global_sub, tmp_subgrid_recon.global_bbox_mmxx);
    tmp_subgrid_skl_global_sub = grid_c_skl_global_sub(tmp_skl_vxl_in_subgrid_Q, :);
    tmp_subgrid_skl_sub = tmp_subgrid_skl_global_sub - tmp_subgrid_recon.global_bbox_mmxx(1:3) + 1;
    tmp_subgrid_skl_ind = sub2ind(tmp_subgrid_recon.block_size, tmp_subgrid_skl_sub(:, 1), ...
        tmp_subgrid_skl_sub(:, 2), tmp_subgrid_skl_sub(:, 3));
    tmp_subgrid_skl_r = grid_c_r(tmp_skl_vxl_in_subgrid_Q);
    % Error: tmp_subgrid_skl_ind is not in the same order as
    % tmp_subgrid_skl_r! - seems have been solved. 
%     assert(all(tmp_subgrid_skl_r == full(vessel_graph.radius(grid_c_skl_ind(tmp_skl_vxl_in_subgrid_Q)))), ...
%         'Mismatch skeleton indices and radius');
    tmp_subgrid_recon.skl_ind = uint32(tmp_subgrid_skl_ind);
    tmp_subgrid_recon.skl_r = tmp_subgrid_skl_r;
    
    % Features
    tmp_subgrid_recon.link.features = vessel_graph.link.features(full(tmp_link_label_2_list_ind(tmp_link_label_in_bbox_cell{iter_int_grid})), :);
    tmp_subgrid_recon.node.features = vessel_graph.node.features(full(tmp_node_label_2_list_ind(tmp_node_label_in_bbox_cell{iter_int_grid})), :);
    tmp_subgrid_recon.link.features.in_bbox_skel_ratio = tmp_link_voxel_in_bbox_ratio_cell{iter_int_grid};
    % In bbox features
    tmp_subgrid_recon.link.features.in_bbox_length = tmp_subgrid_recon.link.features.length .* tmp_subgrid_recon.link.features.in_bbox_skel_ratio;
    tmp_subgrid_recon.link.features.in_bbox_surface_area = tmp_subgrid_recon.link.features.surface_area .* tmp_subgrid_recon.link.features.in_bbox_skel_ratio;
    tmp_subgrid_recon.link.features.in_bbox_volume = tmp_subgrid_recon.link.features.volume .* tmp_subgrid_recon.link.features.in_bbox_skel_ratio;
    
    tmp_subgrid_recon_mask = crop_bbox3(recon_mask, tmp_internal_subgrid_bbox_mmll_grid_c(iter_int_grid, :), 'default');
    tmp_subgrid_recon.ind = uint32(find(tmp_subgrid_recon_mask));
    tmp_subgrid_recon.max_proj_1 = squeeze(max(tmp_subgrid_recon_mask, [], 1));
    tmp_subgrid_recon.max_proj_2 = squeeze(max(tmp_subgrid_recon_mask, [], 2));
    tmp_subgrid_recon.max_proj_3 = max(tmp_subgrid_recon_mask, [], 3);
    
    tmp_subgrid_stat = struct;
    tmp_subgrid_stat.mask_volume = nnz(tmp_subgrid_recon_mask);
    tmp_subgrid_stat.mask_volume_density = tmp_subgrid_stat.mask_volume / prod(tmp_subgrid_recon.block_size);
    tmp_subgrid_recon_mask_erode = imerode(tmp_subgrid_recon_mask, strel('sphere', 1));
    tmp_subgrid_stat.num_surface_voxel = nnz(tmp_subgrid_recon_mask & ~ tmp_subgrid_recon_mask_erode);
    tmp_cc = bwconncomp(tmp_subgrid_recon_mask);
    tmp_cc_surf_area = regionprops3(tmp_cc, 'SurfaceArea');
    % Unknow biased. Maybe over-estimate the surface area because of the edge
    % effect. Maybe still under estimate the surface area due to the
    % finite voxel size.
    tmp_subgrid_stat.mask_surface_area_density_mm2_mm3 = sum(tmp_cc_surf_area.SurfaceArea) ...
        / prod(tmp_subgrid_recon.block_size) * 1e3;
    % Other cube stat
    tmp_subgrid_stat.link_length_density_m_mm3 = sum(tmp_subgrid_recon.link.features.in_bbox_length) / prod(tmp_subgrid_recon.block_size) * 1e3;
    tmp_subgrid_stat.link_surface_area_density_mm2_mm3 = sum(tmp_subgrid_recon.link.features.in_bbox_surface_area) ...
        / prod(tmp_subgrid_recon.block_size) * 1e3;
    tmp_subgrid_stat.link_volume_density = sum(tmp_subgrid_recon.link.features.in_bbox_volume) ...
        / prod(tmp_subgrid_recon.block_size);
    
    tmp_is_capillary_Q = (tmp_subgrid_recon.link.features.dt_median <= opt_comp_graph_feature.capillary_max_radius);
    tmp_subgrid_stat.capillary_length_density_m_mm3 = sum(tmp_subgrid_recon.link.features.in_bbox_length(tmp_is_capillary_Q)) ...
        / prod(tmp_subgrid_recon.block_size) * 1e3;
    tmp_subgrid_stat.capillary_surface_area_density_mm2_mm3 = sum(tmp_subgrid_recon.link.features.in_bbox_surface_area(tmp_is_capillary_Q)) ...
        / prod(tmp_subgrid_recon.block_size) * 1e3;
    tmp_subgrid_stat.capillary_volume_density = sum(tmp_subgrid_recon.link.features.in_bbox_volume(tmp_is_capillary_Q)) ...
        / prod(tmp_subgrid_recon.block_size);
    
    % Capillary vs vessel fraction 
    tmp_subgrid_stat.cap2vsl_length_fraction = tmp_subgrid_stat.capillary_length_density_m_mm3 ./ tmp_subgrid_stat.link_length_density_m_mm3;
    tmp_subgrid_stat.cap2vsl_surf_area_fraction = tmp_subgrid_stat.capillary_surface_area_density_mm2_mm3 ./ tmp_subgrid_stat.link_surface_area_density_mm2_mm3;
    tmp_subgrid_stat.cap2vsl_vol_fraction = tmp_subgrid_stat.capillary_volume_density ./ tmp_subgrid_stat.link_volume_density;   
    
    % Analysis anisotropy
%   Volume-weighted - these two steps might be merged together
    if opt_comp_graph_feature.computeQ.volume_weighted_ep2ep_anisotropy
        if isempty(tmp_subgrid_graph)
            tmp_subgrid_graph = fun_skeleton_to_graph(tmp_subgrid_skl_ind, tmp_subgrid_recon.block_size);
            tmp_subgrid_graph.radius = sparse(double(tmp_subgrid_skl_ind), 1,...
                double(tmp_subgrid_skl_r), prod(tmp_subgrid_recon.block_size), 1);
        end
        tmp_subgrid_stat.anisotropy_all_vw = fun_analysis_get_anisotropy_stat_from_vessel_graph(tmp_subgrid_graph, ...
              [0, inf], 'volume');
          tmp_subgrid_stat.anisotropy_all_lw = fun_analysis_get_anisotropy_stat_from_vessel_graph(tmp_subgrid_graph, ...
              [0, inf], 'length');
    end
    
    if opt_comp_graph_feature.computeQ.capillary_volume_weighted_ep2ep_anisotropy
        if isempty(tmp_subgrid_graph)
            tmp_subgrid_graph = fun_skeleton_to_graph(tmp_subgrid_skl_ind, tmp_subgrid_recon.block_size);
            tmp_subgrid_graph.radius = sparse(double(tmp_subgrid_skl_ind), 1,...
                double(tmp_subgrid_skl_r), prod(tmp_subgrid_recon.block_size), 1);
        end
        tmp_subgrid_stat.anisotropy_capillary_vw = fun_analysis_get_anisotropy_stat_from_vessel_graph(tmp_subgrid_graph, ...
            [0, opt_comp_graph_feature.capillary_max_radius], 'volume');
        tmp_subgrid_stat.anisotropy_capillary_lw = fun_analysis_get_anisotropy_stat_from_vessel_graph(tmp_subgrid_graph, ...
            [0, opt_comp_graph_feature.capillary_max_radius], 'length');
    end
    tmp_subgrid_recon.stat = tmp_subgrid_stat;
    DataManager.write_block_mask(tmp_subgrid_recon, dataset_name, stack, tmp_subgrid_recon.version, ...
        tmp_subgrid_recon.grid_idx_1, tmp_subgrid_recon.grid_idx_2, tmp_subgrid_recon.grid_layer);
end
DataManager.write_graph_in_block(vessel_graph, dataset_name, stack, output_graph_name, ...
    tmp_grid_c_sub(1), tmp_grid_c_sub(2), tmp_grid_c_sub(3));
exit_code = 0;
end