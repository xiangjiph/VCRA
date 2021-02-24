function psf_stat = fun_radius_estimation_get_PSF_stat_in_image_stack(psf_est_int, ...
    tile_image, image_voxel_size_um, skeleton_cc_sub_um, skeleton_cc_r_um)

%% Parameters and parse inputs
psf_est_um = [0.6, 4.0];
ori_num_voxel_half = 2;
ori_num_voxel = ori_num_voxel_half * 2 + 1;

num_link = numel(skeleton_cc_sub_um);
assert(num_link == numel(skeleton_cc_r_um), 'Mismatch array size');

tile_im_size = size(tile_image);
%% Initialization
[link_ori_vec, best_fit_PSF_FWHM] = deal(cell(num_link, 1));
[voxel_dt, refine_r_est, ...
    best_fit_corr] = deal(cell(num_link, 1));
%% Computation
tmp_re_est_tic = tic;
for iter_link = 1 : num_link
    tmp_link_cc_sub_um = skeleton_cc_sub_um{iter_link};
    tmp_link_voxel_r = double(skeleton_cc_r_um{iter_link});
    tmp_num_voxel = numel(tmp_link_voxel_r);
    assert(size(tmp_link_cc_sub_um, 1) == tmp_num_voxel, 'Mismatch array size');
%     tmp_link_r_med = max(1, round(median(tmp_link_voxel_r)));
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
    [tmp_PSF_FWHM, tmp_ori_vec] = deal(nan(3, num_r_est));
    
    for iter_est = 1 : num_r_est
        test_voxel_idx = test_cc_r_est_idx(iter_est);
        % Coordinate of neighbor voxels - for estimating the tilt angle
        tmp_sample_ind_min = max(1, test_voxel_idx - ori_num_voxel_half);
        tmp_sample_ind_max = min(tmp_num_voxel, test_voxel_idx + ori_num_voxel_half);
        tmp_neighbor_voxel_sub_um = tmp_link_cc_sub_um(tmp_sample_ind_min : tmp_sample_ind_max, :);
        % Get the first round of radius estimation
        tmp_dt(iter_est) = tmp_link_voxel_r(test_voxel_idx);        
        % Compute the orientation vector
        tmp_vec = fun_radius_estimation_get_segment_orientation_vector(tmp_neighbor_voxel_sub_um);
        tmp_ori_vec(:, iter_est) = tmp_vec;
        % Find the corresponding point in the full resolution image tile
        tmp_neighbor_voxel_sub_r = round(tmp_neighbor_voxel_sub_um ./ image_voxel_size_um);
        tmp_neighbor_voxel_sub_r = max(1, min(tile_im_size, tmp_neighbor_voxel_sub_r));
                
        est_psf_str = fun_radius_estimation_get_PSF_estimation(psf_est_int, tile_image, ...
            tmp_neighbor_voxel_sub_r, tmp_dt(iter_est), tmp_vec, ...
            image_voxel_size_um, psf_est_um, false);
        
        if ~isempty(est_psf_str)
            tmp_est_r(iter_est) = est_psf_str.vsl_r_est;
            tmp_PSF_FWHM(:, iter_est) = est_psf_str.best_fit_FWHM';
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
    link_ori_vec{iter_link} = tmp_ori_vec.';
    voxel_dt{iter_link} = tmp_dt;
    refine_r_est{iter_link} = tmp_est_r;
    best_fit_PSF_FWHM{iter_link} = tmp_PSF_FWHM.';
    best_fit_corr{iter_link} = tmp_fit_corr;
    fprintf('Finish estimating PSF near link %d\n', iter_link);
end
fprintf('Finish estimating the PSF. Elapsed time is %f seconds.\n', toc(tmp_re_est_tic));
%% Return data
psf_stat = struct;
psf_stat.link_ori_vec = link_ori_vec;
psf_stat.best_fit_PSF_FWHM = best_fit_PSF_FWHM;
psf_stat.voxel_dt = voxel_dt;
psf_stat.refine_r_est = refine_r_est;
psf_stat.best_fit_corr = best_fit_corr;
end