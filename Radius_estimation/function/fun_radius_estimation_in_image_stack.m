function [est_radius] = fun_radius_estimation_in_image_stack(est_psf_int, tile_image, image_voxel_size_um, ...
    skeleton_cc_sub_um, skeleton_cc_r_um)
% This function is for post-perfusion radius estimation (refinement) 

%% Parameters
% The maximum size of vessel to estimate the radius
max_r_to_est_um = max(est_psf_int.radius(:));
ori_num_voxel_half = 2;
ori_num_voxel = ori_num_voxel_half * 2 + 1;
% Number of iteration for radius estimation
num_r_est_round = 6;
% If the estimated radius is 2 um larger than the original estimation,
% reject the update and use the origninal estimation
max_est_update_diff = 2;

% Debug setting
vis_Q = false;
%% Parse input 
num_cc = numel(skeleton_cc_sub_um);
assert(num_cc == numel(skeleton_cc_r_um), 'Mismatch input cc array size');
tile_image_size = size(tile_image);
%% Initialization
[est_link_is_valid_Q, est_link_radius, est_link_ori_vec_z] = deal(cell(num_cc, 1));
est_radius_stat = struct;
[est_radius_stat.median, est_radius_stat.mean, est_radius_stat.std] = deal(nan(num_cc, 1));
%% Computation
tmp_cube_tic = tic;
for iter_cc = 1 : num_cc
    tmp_cc_sub_um = skeleton_cc_sub_um{iter_cc};
    tmp_cc_r = double(skeleton_cc_r_um{iter_cc});
    tmp_cc_num_voxel = size(tmp_cc_sub_um, 1);
    assert(tmp_cc_num_voxel == numel(tmp_cc_r));
    
    tmp_cc_r_med = median(tmp_cc_r);
%     tmp_cc_half_step = ceil(tmp_cc_r_med / 2 );
    tmp_cc_half_step = 0;
   
    if tmp_cc_r_med <= max_r_to_est_um
        % Determine the broadcasting indices 
        if tmp_cc_num_voxel >= ori_num_voxel
            test_cc_r_est_idx = (ori_num_voxel_half + 1) : (2 * tmp_cc_half_step + 1) : (tmp_cc_num_voxel - ori_num_voxel_half);
            test_cc_r_seg_upper_ind = test_cc_r_est_idx + tmp_cc_half_step;
            test_cc_r_seg_upper_ind(end) = tmp_cc_num_voxel;
            test_cc_r_seg_num_vxl = diff([0, test_cc_r_seg_upper_ind]);
            assert(sum(test_cc_r_seg_num_vxl) == tmp_cc_num_voxel, 'The total length of the segments does not equal the length of the vessel');
        else
            test_cc_r_est_idx = ceil(tmp_cc_num_voxel / 2);
            test_cc_r_seg_num_vxl = tmp_cc_num_voxel;
        end
        
        num_r_est = numel(test_cc_r_est_idx);
        tmp_cc_est_r_list = nan(num_r_est, 1);
        tmp_cc_ori_vec_z = nan(num_r_est, 1);
        for iter_est = 1 : num_r_est
            tmp_voxel_idx = test_cc_r_est_idx(iter_est);  
            tmp_voxel_r = tmp_cc_r(tmp_voxel_idx);
            % Coordinate of neighbor voxels - for estimating the tilt angle
            tmp_sample_ind_min = max(1, tmp_voxel_idx - ori_num_voxel_half);
            tmp_sample_ind_max = min(tmp_cc_num_voxel, tmp_voxel_idx + ori_num_voxel_half);
            tmp_neighbor_voxel_sub_um = tmp_cc_sub_um(tmp_sample_ind_min : tmp_sample_ind_max, :);
            % Estimate the segmentation orientation
            tmp_seg_ori_vec = fun_radius_estimation_get_segment_orientation_vector(tmp_neighbor_voxel_sub_um);
            % Find the corresponding point in the full resolution image tile
            tmp_neighbor_voxel_sub_r = round(tmp_neighbor_voxel_sub_um ./ image_voxel_size_um);
            tmp_neighbor_voxel_sub_r = max(1, min(tile_image_size, tmp_neighbor_voxel_sub_r));
            % Visualization
            if vis_Q && false
                tmp_voxel_sub = tmp_cc_sub_um(tmp_voxel_idx, :);
                test_voxel_sub_r = round(tmp_voxel_sub ./ image_voxel_size_um);
                fig_hdl = figure;
                ax_2 = axes(fig_hdl);
                imagesc(ax_2, tile_image(:, :, test_voxel_sub_r(3)));
                hold(ax_2, 'on');
                scatter(ax_2, test_voxel_sub_r(2), test_voxel_sub_r(1), 25, 'LineWidth', 4);
                ax_2.DataAspectRatio = [1,1,1];
                ax_2.Visible = 'on';
                ax_2.Title.String = 'Rendered image';
                ax_2.XAxis.Visible = 'off';
                ax_2.YAxis.Visible = 'off';
                colormap('jet');
                colorbar(ax_2);
            end
            
            tmp_r_est = fun_radius_estimation_by_adp_thrld_DT(tile_image, tmp_neighbor_voxel_sub_r,...
                tmp_voxel_r, tmp_seg_ori_vec, ...
                est_psf_int.n_min_edge_int, image_voxel_size_um, num_r_est_round);
            
            if ~isnan(tmp_r_est)
                tmp_r_est_correction = abs(tmp_r_est - tmp_voxel_r);
                if tmp_r_est_correction < max_est_update_diff
                    tmp_cc_est_r_list(iter_est) = tmp_r_est;
                end
                tmp_cc_ori_vec_z(iter_est) = abs(tmp_seg_ori_vec(3));
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
                    tmp_cc_ori_vec_z(tmp_is_nan_Q) = median(tmp_cc_ori_vec_z(~tmp_is_nan_Q));
                end
            end
            % Use the estimated radius in selected skeleton voxels as the
            % radius of its neighboring skeleton voxels - duplicating
            tmp_cc_est_r_list = repelem(tmp_cc_est_r_list, test_cc_r_seg_num_vxl, 1);
            tmp_cc_ori_vec_z = repelem(tmp_cc_ori_vec_z, test_cc_r_seg_num_vxl, 1);
        else
            % If all the estimated radius is nan, use the original
            % estimation
            tmp_cc_est_r_list = tmp_cc_r;
            tmp_cc_ori_vec_z = nan(size(tmp_cc_r));
        end              
        tmp_is_nan_Q = repelem(tmp_is_nan_Q, test_cc_r_seg_num_vxl, 1);
    else
        tmp_cc_est_r_list = tmp_cc_r;
        tmp_is_nan_Q = false(size(tmp_cc_r));
        tmp_cc_ori_vec_z = nan(size(tmp_cc_r));
    end
    %% Save cc estimation result
    est_link_is_valid_Q{iter_cc} = ~tmp_is_nan_Q;
    % Evaluate the reliability of radius estimation before saving
    est_link_ori_vec_z{iter_cc} = tmp_cc_ori_vec_z;
    est_link_radius{iter_cc} = tmp_cc_est_r_list;
    est_radius_stat.median(iter_cc) = median(tmp_cc_est_r_list, 'omitnan');
    est_radius_stat.mean(iter_cc) = mean(tmp_cc_est_r_list, 'omitnan'); 
    est_radius_stat.std(iter_cc) = std(tmp_cc_est_r_list, 'omitnan');
end
fprintf('Finish estimating radius for the entire block. Elapsed time is %f seconds\n', toc(tmp_cube_tic));
%% Output
est_radius = struct;
est_radius.r = est_link_radius;
est_radius.local_ori_vec_z = est_link_ori_vec_z;
est_radius.valid_Q = est_link_is_valid_Q;
est_radius.cc_stat = est_radius_stat;
end