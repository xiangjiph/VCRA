function exit_code = fun_radius_estimation_in_combined_grid(re_opt, grid_c_list_ind)


%% Initialization
dataset_name = re_opt.dataset_name;
stack = re_opt.stack;
load_skel_version = re_opt.load_skel_version;
write_skel_version = re_opt.write_skel_version;

persistent DataManager;
if isempty(DataManager)
    DataManager = FileManager;
end

grid_info = re_opt.grid_info;
data_octree = grid_info.octree;
grid_voxel_size = grid_info.voxel_size_um;
% Running setting
parallel_load_Q = re_opt.parallel_load_Q;

% Radius parameter setting
re_grid = re_opt.grid_radius_estimation;
est_psf_int = re_opt.psf_edge_int_interpolation;
% The maximum size of vessel to estimate the radius
max_r_to_est_um = re_opt.max_r_to_est_um; 
ori_num_voxel_half = re_opt.num_ori_vec_vxl_half;
% Number of iteration for radius estimation
num_r_est_round = re_opt.num_iteration;
% If the estimated radius is 2 um larger than the original estimation,
% reject the update and use the origninal estimation
max_est_update_diff = re_opt.max_ext_to_ori_r_um;

% Debug setting
vis_Q = false;
%% Load skeleton files to generate the graph
data_info = fun_graph_get_grid_c_block_info(re_grid, grid_c_list_ind);

grid_c_array_ind = re_grid.bbox_grid_ind_list(grid_c_list_ind, :);
grid_c_subgrid_array_sub = re_grid.sub_grid_sub{grid_c_array_ind};
grid_c_subgrid_global_bbox_mmxx_pxl = re_grid.sub_grid_bbox_mmxx{grid_c_array_ind};
num_subgrid = size(grid_c_subgrid_array_sub, 1);

grid_c_bbox_xyz_pxl_mmxx = re_grid.bbox_xyz_mmxx_pixel_list(grid_c_list_ind, :);
grid_c_bbox_xyz_pxl_ll = grid_c_bbox_xyz_pxl_mmxx(4:6) - grid_c_bbox_xyz_pxl_mmxx(1:3) + 1;

grid_c_skl_ind = cell(num_subgrid, 1);
grid_c_skl_r = cell(num_subgrid, 1);

tmp_tic = tic;
% subgrid_est_bg = nan(num_subgrid, 1);
% subgrid_est_bg_std = nan(num_subgrid, 1);
for iter_subgrid = 1 : num_subgrid
    tmp_idx_1 = grid_c_subgrid_array_sub(iter_subgrid, 1);
    tmp_idx_2 = grid_c_subgrid_array_sub(iter_subgrid, 2);
    tmp_layer = grid_c_subgrid_array_sub(iter_subgrid, 3);
    try
        tmp_skel_str = DataManager.load_block_skl(dataset_name, stack, ...
            load_skel_version, tmp_idx_1, tmp_idx_2, tmp_layer);
% Might be useful for background estimation
%         tmp_seg_info = load(DataManager.fp_block_mask_file(dataset_name, stack, ...
%             '240_cube', tmp_idx_1, tmp_idx_2, tmp_layer), 'record');
%         subgrid_est_bg(iter_subgrid) = tmp_seg_info.record.raw.est_bg;
%         subgrid_est_bg_std(iter_subgrid) = tmp_seg_info.record.raw.est_bg_std;
    catch ME
        fprintf('Fail to load the skeleton file of cube (%d, %d, %d). Skip this block\n', tmp_idx_1, tmp_idx_2, tmp_layer);
        fprintf('Error message is:\n%s', ME.message);
        continue;
    end
    % Conver to combined grid coordinate
    tmp_skl_local_sub = fun_ind2sub(tmp_skel_str.block_size, tmp_skel_str.ind);
    tmp_skl_global_sub = tmp_skl_local_sub + tmp_skel_str.global_bbox_mmll(1:3) - 1;
    tmp_skl_cg_sub = tmp_skl_global_sub - grid_c_bbox_xyz_pxl_mmxx(1:3) + 1;
    tmp_skl_cg_ind = sub2ind(grid_c_bbox_xyz_pxl_ll, tmp_skl_cg_sub(:, 1), ...
        tmp_skl_cg_sub(:, 2), tmp_skl_cg_sub(:, 3));
    % Save to the cell array
    grid_c_skl_ind{iter_subgrid} = tmp_skl_cg_ind;
    grid_c_skl_r{iter_subgrid} = double(tmp_skel_str.r);
end
fprintf('Finish reading skeleton data. Elapsed time is %f seconds.\n', toc(tmp_tic));
% Merge cell array
grid_c_skl_ind = cat(1, grid_c_skl_ind{:});
grid_c_skl_r = cat(1, grid_c_skl_r{:});
assert(numel(grid_c_skl_ind) == numel(grid_c_skl_r), 'Number of skeleton voxel does not equal the number of radius data');
% Remove duplicated voxels in the overlapping region
[grid_c_skl_ind, tmp_unique_ind, ~] = unique(grid_c_skl_ind, 'stable');
grid_c_skl_r = grid_c_skl_r(tmp_unique_ind);

num_skl_voxel = numel(grid_c_skl_ind);
if num_skl_voxel == 0
    fprintf('Empty block. Return.\n');
    exit_code = 4;
    return;
end
% Convert to graph
vessel_graph = fun_skeleton_to_graph(grid_c_skl_ind, grid_c_bbox_xyz_pxl_ll);
vessel_graph.radius = sparse(grid_c_skl_ind, 1, grid_c_skl_r, prod(grid_c_bbox_xyz_pxl_ll), 1);
vessel_graph.info = data_info;
fprintf('Finish generating vessel graph\n');
%% Load full resolution images
% Determine the bounding box for the rendered data
% Expand the bounding box by 1 voxel on xy direction to avoid rounding error
load_im_ds_bbox_mmxx = grid_c_bbox_xyz_pxl_mmxx;
load_im_ds_bbox_mmxx(4:5) = min(re_grid.data_size(1:2), grid_c_bbox_xyz_pxl_mmxx(4:5) + 1);

tmp_tic = tic;
[tile_image, tile_image_bbox_r_mm] = fun_radius_estimation_load_rendered_image(...
    data_octree, load_im_ds_bbox_mmxx, grid_voxel_size, parallel_load_Q);
fprintf('Finish loading rendered image. Elapsed time is %f seconds\n', toc(tmp_tic));
if isfield(re_opt, 'render_data_medfilt3_Q') && re_opt.render_data_medfilt3_Q
    tmp_tic = tic;
    tile_image = medfilt3(tile_image);
    fprintf('Finish applying 3D median filter. Elapsed time is %f seconds\n', toc(tmp_tic));
end
%% Estimate the radius for all the link segment in 240-cube
[est_link_is_valid_Q, est_link_radius] = deal(cell(vessel_graph.link.num_cc, 1));

est_radius_stat = struct;
[est_radius_stat.median, est_radius_stat.mean, est_radius_stat.std] = deal(nan(vessel_graph.link.num_cc, 1));
ori_num_voxel = ori_num_voxel_half * 2 + 1;

tmp_cube_tic = tic;
for test_link_label = 1 : vessel_graph.link.num_cc
    test_cc_ind = vessel_graph.link.cc_ind{test_link_label};
    test_cc_r = full(vessel_graph.radius(test_cc_ind));
    test_cc_r_med = median(test_cc_r);
    test_cc_half_step = ceil(test_cc_r_med / 2 );
    
    if test_cc_r_med <= max_r_to_est_um
        test_cc_num_voxel = numel(test_cc_ind);
        test_cc_sub = fun_ind2sub(vessel_graph.num.mask_size, test_cc_ind);
        
        if test_cc_num_voxel >= ori_num_voxel
            test_cc_r_est_idx = (ori_num_voxel_half + 1) : (2 * test_cc_half_step + 1) : (test_cc_num_voxel - ori_num_voxel_half);
            test_cc_r_seg_upper_ind = test_cc_r_est_idx + test_cc_half_step;
            test_cc_r_seg_upper_ind(end) = test_cc_num_voxel;
            test_cc_r_seg_num_vxl = diff([0, test_cc_r_seg_upper_ind]);
            assert(sum(test_cc_r_seg_num_vxl) == test_cc_num_voxel, 'The total length of the segments does not equal the length of the vessel');
        else
            test_cc_r_est_idx = ceil(test_cc_num_voxel / 2);
            test_cc_r_seg_num_vxl = test_cc_num_voxel ;
        end
        num_r_est = numel(test_cc_r_est_idx);
        tmp_cc_est_r_list = nan(num_r_est, 1);
        for iter_est = 1 : num_r_est
            test_voxel_idx = test_cc_r_est_idx(iter_est);
            % Coordinate of the voxel            
            test_voxel_sub = test_cc_sub(test_voxel_idx, :);
            % Coordinate of neighbor voxels - for estimating the tilt angle
            tmp_sample_ind_min = max(1, test_voxel_idx - ori_num_voxel_half);
            tmp_sample_ind_max = min(test_cc_num_voxel, test_voxel_idx + ori_num_voxel_half);
            
            tmp_neighbor_voxel_sub = test_cc_sub(tmp_sample_ind_min : tmp_sample_ind_max, :);
            % Estimate the segmentation orientation
            tmp_seg_ori_vec = fun_radius_estimation_get_segment_orientation_vector(tmp_neighbor_voxel_sub);
            % Find the corresponding point in the full resolution image tile
            tmp_neighbor_voxel_sub_r = fun_radius_estimation_voxel_local_sub_ds2r(tmp_neighbor_voxel_sub, ...
                grid_c_bbox_xyz_pxl_mmxx(1:3), grid_voxel_size, data_octree.voxel_size, tile_image_bbox_r_mm);
            % Deal with corner case: when tmp_neighbor_voxel_sub_r is
            % outside the bounding box - can happen at the end of the
            % sample 
            tmp_neighbor_voxel_sub_r = max([1,1,1], min(tmp_neighbor_voxel_sub_r, size(tile_image)));
            % Visualization
            if vis_Q && false
                fig_hdl = figure;
                fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
                ax_1 = subplot(1,2,1);
                imagesc(ax_1, grid_c_im(:, :, test_voxel_sub(3)));
                hold(ax_1, 'on');
                ax_1.DataAspectRatio = [1,1,1];
                ax_1.XAxis.Visible = 'off';
                ax_1.YAxis.Visible = 'off';
                ax_1.Title.String = 'Downsampled image';
                scatter(ax_1, test_voxel_sub(2), test_voxel_sub(1), 100, 'LineWidth', 4);
                ax_2 = subplot(1,2,2);
                imagesc(ax_2, tile_image(:, :, test_voxel_sub_tile_im(3)));
                hold(ax_2, 'on');
                scatter(ax_2, test_voxel_sub_tile_im(2), test_voxel_sub_tile_im(1), 100, 'LineWidth', 4);
                ax_2.DataAspectRatio = [1,1,1];
                ax_2.Visible = 'on';
                ax_2.Title.String = 'Rendered image';
                ax_2.XAxis.Visible = 'off';
                ax_2.YAxis.Visible = 'off';
                colormap('gray');
            end
            
            tmp_r_est = fun_radius_estimation_by_adp_thrld_DT(tile_image, tmp_neighbor_voxel_sub_r,...
                test_cc_r(test_voxel_idx), tmp_seg_ori_vec, ...
                est_psf_int.n_min_edge_int, data_octree.voxel_size, num_r_est_round);
            
            if ~isnan(tmp_r_est)
                tmp_r_est_correction = tmp_r_est - test_cc_r(test_voxel_idx);
                if tmp_r_est_correction < max_est_update_diff
                    tmp_cc_est_r_list(iter_est) = tmp_r_est;
                end
            end
        end
       %% Duplicate the radius estimations
        tmp_is_nan_Q = isnan(tmp_cc_est_r_list);
        if ~all(tmp_is_nan_Q)
            if any(tmp_is_nan_Q)
                % If part of the estimated radius are nan, replace them
                % with the median raidus. 
                tmp_est_r_not_nan = tmp_cc_est_r_list(~tmp_is_nan_Q);
                if ~isempty(tmp_est_r_not_nan)
                    tmp_cc_est_r_list(tmp_is_nan_Q) = median(tmp_cc_est_r_list(~tmp_is_nan_Q));
                end
            end
            % Use the estimated radius in selected skeleton voxels as the
            % radius of its neighboring skeleton voxels - duplicating
            tmp_cc_est_r_list = repelem(tmp_cc_est_r_list, test_cc_r_seg_num_vxl, 1);
        else
            % If all the estimated radius is nan, use the original
            % estimation
            tmp_cc_est_r_list = test_cc_r;
        end              
        tmp_is_nan_Q = repelem(tmp_is_nan_Q, test_cc_r_seg_num_vxl, 1);
    else
        tmp_cc_est_r_list = test_cc_r;
        tmp_is_nan_Q = false(size(test_cc_r));
    end
    est_link_is_valid_Q{test_link_label} = ~tmp_is_nan_Q;
    % Evaluate the reliability of radius estimation before saving
    est_link_radius{test_link_label} = tmp_cc_est_r_list;
    est_radius_stat.median(test_link_label) = median(tmp_cc_est_r_list, 'omitnan');
    est_radius_stat.mean(test_link_label) = mean(tmp_cc_est_r_list, 'omitnan'); 
    est_radius_stat.std(test_link_label) = std(tmp_cc_est_r_list, 'omitnan');
end
fprintf('Finish estimating radius for the entire block. Elapsed time is %f seconds\n', toc(tmp_cube_tic));
%% Get the radius for the node voxels
node_r = cell(vessel_graph.node.num_cc, 1);
for iter_node = 1 : vessel_graph.node.num_cc
    tmp_connected_link_label = vessel_graph.node.connected_link_label{iter_node};
    tmp_connected_link_med_r = median(est_radius_stat.median(tmp_connected_link_label));
    tmp_connected_link_med_r = repelem(tmp_connected_link_med_r, vessel_graph.node.num_voxel_per_cc(iter_node), 1);
    node_r{iter_node} = tmp_connected_link_med_r;
end
node_r = cat(1, node_r{:});
node_valid_Q = true(size(node_r));
%% Update radius and write to disk
update_ind = cat(1, vessel_graph.link.pos_ind, vessel_graph.node.pos_ind);
% Deal with the rest of the voxels
non_node_link_voxels = cat(1, vessel_graph.isoloop.pos_ind, vessel_graph.isopoint.pos_ind);
non_node_link_r = full(vessel_graph.radius(non_node_link_voxels));
non_node_link_valid_Q = false(size(non_node_link_r));

update_ind = cat(1, update_ind, non_node_link_voxels);
update_r = cat(1, est_link_radius{:}, node_r, non_node_link_r);
update_valid_Q = cat(1, est_link_is_valid_Q{:}, node_valid_Q, non_node_link_valid_Q);
assert(numel(update_ind) == numel(update_r), 'Number of element mismatches');
assert(numel(update_ind) == numel(grid_c_skl_ind), 'Number of updated voxels does not match the number of voxles in the graph');
assert(numel(update_valid_Q) == numel(update_ind), 'Number of updated voxels does not match the number of element in the valid_Q vector');

update_sub = fun_ind2sub(vessel_graph.num.mask_size, update_ind);
update_sub_global = bsxfun(@plus, update_sub, grid_c_bbox_xyz_pxl_mmxx(1:3) - 1);

update_skel_str_0 = struct;
update_skel_str_0.dataset_name = dataset_name;
update_skel_str_0.stack = stack;
update_skel_str_0.grid_name = write_skel_version;
update_skel_str_0.dataset_size = grid_info.data_size;
[update_skel_str_0.global_bbox_mmxx, update_skel_str_0.global_bbox_mmll, ...
    update_skel_str_0.idx_1, update_skel_str_0.idx_2, update_skel_str_0.layer, ...
    update_skel_str_0.block_size, update_skel_str_0.ind, update_skel_str_0.r] = deal([]);

for iter_subgrid = 1 : num_subgrid
    tmp_idx_1 = grid_c_subgrid_array_sub(iter_subgrid, 1);
    tmp_idx_2 = grid_c_subgrid_array_sub(iter_subgrid, 2);
    tmp_layer = grid_c_subgrid_array_sub(iter_subgrid, 3);   
    
    tmp_update_skel_str = update_skel_str_0;
    tmp_update_skel_str.idx_1 = tmp_idx_1;
    tmp_update_skel_str.idx_2 = tmp_idx_2;
    tmp_update_skel_str.layer = tmp_layer;    
    tmp_update_skel_str.global_bbox_mmxx = grid_c_subgrid_global_bbox_mmxx_pxl(iter_subgrid, :);
    tmp_update_skel_str.block_size = tmp_update_skel_str.global_bbox_mmxx(4:6) - tmp_update_skel_str.global_bbox_mmxx(1:3) + 1;
    tmp_update_skel_str.global_bbox_mmll = tmp_update_skel_str.global_bbox_mmxx;
    tmp_update_skel_str.global_bbox_mmll(4:6) = tmp_update_skel_str.block_size;
    % Select the voxels inside the global bbox mmxx
    tmp_is_in_bbox_Q = fun_voxel_sub_in_bbox_mmxx_Q(update_sub_global, tmp_update_skel_str.global_bbox_mmxx);
    tmp_local_sub = update_sub_global(tmp_is_in_bbox_Q, :) - tmp_update_skel_str.global_bbox_mmxx(1:3) + 1;
    tmp_local_ind = sub2ind(tmp_update_skel_str.block_size, tmp_local_sub(:, 1), tmp_local_sub(:, 2), tmp_local_sub(:, 3));
    
    tmp_update_skel_str.ind = uint32(tmp_local_ind);
    tmp_update_skel_str.r = single(update_r(tmp_is_in_bbox_Q));
    tmp_update_skel_str.valid_Q = update_valid_Q(tmp_is_in_bbox_Q);
    % Write data
    DataManager.write_block_skl_file(tmp_update_skel_str, tmp_update_skel_str.dataset_name, ...
        tmp_update_skel_str.stack, tmp_update_skel_str.grid_name, ...
        tmp_update_skel_str.idx_1, tmp_update_skel_str.idx_2, tmp_update_skel_str.layer);
end
fprintf('Finish writing skeleton structure with updated radius estimation to hard drive\n');
%% Save graph
vessel_graph.radius_refined = sparse(update_ind, 1, update_r, prod(vessel_graph.num.mask_size), 1);
vessel_graph.link.features.r_est_med = est_radius_stat.median;
vessel_graph.link.features.r_est_mean = est_radius_stat.mean;
vessel_graph.link.features.r_est_std = est_radius_stat.std;
DataManager.write_graph_in_block(vessel_graph, dataset_name, stack, vessel_graph.info.grid_version, ...
    vessel_graph.info.idx_1, vessel_graph.info.idx_2, vessel_graph.info.layer);
exit_code = 0;
end