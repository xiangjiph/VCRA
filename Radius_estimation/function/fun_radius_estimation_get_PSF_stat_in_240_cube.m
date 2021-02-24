function psf_stat = fun_radius_estimation_get_PSF_stat_in_240_cube(psf_est_int, grid_info, skel_version, cube_grid_ind, medfilt_raw_Q)

if nargin < 5
    medfilt_raw_Q = false;
end

persistent DataManager
if isempty(DataManager)
    DataManager = FileManager;
end

dataset_name = grid_info.dataset_name;
stack = grid_info.stack;
data_octree = grid_info.octree;
grid_voxel_size = grid_info.voxel_size_um;
%% Load 240 cube
cube_grid_sub = grid_info.bbox_grid_sub_list(cube_grid_ind, :);

skl_str = DataManager.load_block_skl(dataset_name, stack, skel_version, ...
    cube_grid_sub(1), cube_grid_sub(2), cube_grid_sub(3));
cube_global_mmxx_pxl = skl_str.global_bbox_mmxx;


vessel_graph = fun_skeleton_to_graph(skl_str.ind, skl_str.block_size);
vessel_graph.radius = sparse(double(skl_str.ind), 1, double(skl_str.r), prod(skl_str.block_size), 1);

vessel_graph.link.features.dt_median = nan(vessel_graph.link.num_cc, 1);
for iter_cc = 1 : vessel_graph.link.num_cc
    tmp_ind = vessel_graph.link.cc_ind{iter_cc};
    vessel_graph.link.features.dt_median(iter_cc) = median(full(vessel_graph.radius (tmp_ind)));
end
%% Debug
% vessel_im_ds = DataManager.load_block_data(dataset_name, stack, '240_cube', ...
%     cube_grid_sub(1), cube_grid_sub(2), cube_grid_sub(3));
%% Use the bounding box of the image to determine the bounding box of the
% rendered data
% Expand the bounding box
skel_bbox_mmxx = skl_str.global_bbox_mmxx;
exp_length = 10;
skel_bbox_mmxx(1:2) = max(1, skel_bbox_mmxx(1:2) - exp_length);
skel_bbox_mmxx(4:5) = min(grid_info.data_size(1:2), skel_bbox_mmxx(4:5) + exp_length);

tmp_tic = tic;
[tile_image, tile_image_bbox_r_mm] = fun_radius_estimation_load_rendered_image(...
    data_octree, skel_bbox_mmxx, grid_voxel_size, false);
fprintf('Finish loading rendered image. Elapsed time is %f seconds\n', toc(tmp_tic));
if medfilt_raw_Q 
    tile_image = medfilt3(tile_image);
    fprintf('Finish applying 3D median filter\n');
end
%% PSF - radius joint estimation
ori_num_voxel_half = 2;
ori_num_voxel = ori_num_voxel_half * 2 + 1;

num_link = vessel_graph.link.num_cc;
[link_ori_vec, best_fit_PSF_FWHM] = deal(cell(num_link, 1));
[voxel_dt, refine_r_est, ...
    best_fit_corr] = deal(cell(num_link, 1));
tmp_re_est_tic = tic;
for iter_link = 1 : num_link
    tmp_link_label = iter_link;
    
    tmp_link_cc_ind = vessel_graph.link.cc_ind{tmp_link_label};
    tmp_num_voxel = numel(tmp_link_cc_ind);
    tmp_link_cc_sub = fun_ind2sub(vessel_graph.num.mask_size, tmp_link_cc_ind);
    tmp_link_voxel_r = double(full(vessel_graph.radius(tmp_link_cc_ind)));
%     tmp_link_r_med = ceil(median(tmp_link_voxel_r));
%     tmp_half_step = ceil(tmp_link_r_med / 2);
%     tmp_est_step = 2 * tmp_half_step + 1;
    tmp_est_step = 2;
    
    if tmp_num_voxel >= ori_num_voxel
        test_cc_r_est_idx = (ori_num_voxel_half + 1) : (tmp_est_step) : (tmp_num_voxel - ori_num_voxel_half);
    else
        test_cc_r_est_idx = ceil(tmp_num_voxel / 2);
    end
    num_r_est = numel(test_cc_r_est_idx);
    
    [tmp_dt, tmp_est_r, tmp_fit_corr] = deal(nan(num_r_est, 1));
    [tmp_PSF_FWHM, tmp_ori_vec] = deal(nan(num_r_est, 3));
    
    for iter_est = 1 : num_r_est
        test_voxel_idx = test_cc_r_est_idx(iter_est);
        % Coordinate of neighbor voxels - for estimating the tilt angle
        tmp_sample_ind_min = max(1, test_voxel_idx - ori_num_voxel_half);
        tmp_sample_ind_max = min(tmp_num_voxel, test_voxel_idx + ori_num_voxel_half);
        tmp_neighbor_voxel_sub = tmp_link_cc_sub(tmp_sample_ind_min : tmp_sample_ind_max, :);
        % Get the first round of radius estimation
        tmp_dt(iter_est) = tmp_link_voxel_r(test_voxel_idx);        
        % Compute the orientation vector
        tmp_vec = fun_radius_estimation_get_segment_orientation_vector(tmp_neighbor_voxel_sub);
        tmp_ori_vec(iter_est, :) = tmp_vec;
        % Find the corresponding point in the full resolution image tile
        tmp_neighbor_voxel_sub_r = fun_radius_estimation_voxel_local_sub_ds2r(tmp_neighbor_voxel_sub, ...
            cube_global_mmxx_pxl(1:3), grid_voxel_size, data_octree.voxel_size, tile_image_bbox_r_mm);
        if abs(tmp_vec(3)) <= sind(60) && tmp_dt(iter_est) <= 3 && num_r_est > 20 && ...
                iter_est == 11
            fprintf('debug\n');
        end
        est_psf_str = fun_radius_estimation_get_PSF_estimation(psf_est_int, tile_image, ...
            tmp_neighbor_voxel_sub_r, tmp_dt(iter_est), tmp_vec, ...
            data_octree.voxel_size);
        
        if ~isempty(est_psf_str)
            tmp_est_r(iter_est) = est_psf_str.vsl_r_est;
            tmp_PSF_FWHM(iter_est, :) = est_psf_str.best_fit_FWHM';
            tmp_fit_corr(iter_est) = est_psf_str.corr_max;
        end
        %% Debug
%         debug_ds_sub = tmp_neighbor_voxel_sub(3, :);
%         debug_render_sub = tmp_neighbor_voxel_sub_r(3, :);
%         fig_hdl = figure;
%         fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2, 1];
%         ax_hdl_1 = subplot(1,2,1);
%         imagesc(ax_hdl_1, vessel_im_ds(:, :, debug_ds_sub(3)));
%         hold(ax_hdl_1, 'on');
%         scatter(ax_hdl_1, debug_ds_sub(2), debug_ds_sub(1), 'r*');
%         ax_hdl_1.Children(1).SizeData = 100;
%         ax_hdl_1.DataAspectRatio = [1,1,1];
%         ax_hdl_2 = subplot(1,2,2);
%         imagesc(ax_hdl_2, tile_image(:, :, debug_render_sub(3)));
%         hold(ax_hdl_2, 'on');
%         scatter(ax_hdl_2, debug_render_sub(2), debug_render_sub(1), 'r*');
%         ax_hdl_2.Children(1).SizeData = 100;
%         ax_hdl_2.DataAspectRatio = [1,1,1];
    end
    % Save
    link_ori_vec{iter_link} = tmp_ori_vec;
    voxel_dt{iter_link} = tmp_dt;
    refine_r_est{iter_link} = tmp_est_r;
    best_fit_PSF_FWHM{iter_link} = tmp_PSF_FWHM;
    best_fit_corr{iter_link} = tmp_fit_corr;
    fprintf('Finish estimating PSF near link %d\n', iter_link);
end
fprintf('Finish estimating the PSF. Elapsed time is %f seconds.\n', toc(tmp_re_est_tic));
%% Return data
psf_stat = struct;
psf_stat.link_ori_vec = cat(1, link_ori_vec{:});
psf_stat.best_fit_PSF_FWHM = cat(1, best_fit_PSF_FWHM{:});
psf_stat.voxel_dt = cat(1, voxel_dt{:});
psf_stat.refine_r_est = cat(1, refine_r_est{:});
psf_stat.best_fit_corr = cat(1, best_fit_corr{:});
end